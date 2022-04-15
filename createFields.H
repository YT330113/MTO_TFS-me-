Info<< "Reading optimization parameters\n" << endl;
IOdictionary optProperties
(
   IOobject
   (
       "optProperties",
       runTime.constant(),
       mesh,
       IOobject::MUST_READ_IF_MODIFIED,
       IOobject::NO_WRITE
   )
);
/*优化参数*/
scalar raa0(readScalar(optProperties.lookup("raa0")));
scalar mma_init(readScalar(optProperties.lookup("mma_init")));
scalar mma_dec(readScalar(optProperties.lookup("mma_dec")));
scalar mma_inc(readScalar(optProperties.lookup("mma_inc")));
scalar movlim(readScalar(optProperties.lookup("movlim")));

scalar voluse(readScalar(optProperties.lookup("voluse")));  // 体积分数
scalar GeoDim(readScalar(optProperties.lookup("GeoDim"))); // 模型维度，三维/二维

scalar PowerDiss0(readScalar(optProperties.lookup("PowerDiss0"))); // 能量耗散归一化值
scalar PowerDissMax(readScalar(optProperties.lookup("PowerDissMax"))); // 最大能量耗散倍数
scalar PowerDissRelax(readScalar(optProperties.lookup("PowerDissRelax"))); // 能量耗散放倍数

scalar C0(readScalar(optProperties.lookup("C0")));// 结构柔度归一化值
scalar CMax(readScalar(optProperties.lookup("CMax")));// 最大结构柔度倍数
scalar CRelax(readScalar(optProperties.lookup("CRelax")));// 结构柔度放松倍数

scalar filterR(readScalar(optProperties.lookup("filterR"))); //密度过滤半径
dimensionedScalar b("b", dimensionSet(0,-2,0,0,0,0,0),1.0);//filter radius of the PDE filter，初始化PDE过滤半径，在“PDE filter”和“PDE 过滤的链式求导”里用到
dimensionedScalar b1("b1", dimensionSet(0,2,-2,-2,0,0,0),1.0); //在求Ub和fsenshMeanT的时候用到

Switch solid_area(optProperties.lookup("solid_area"));
Switch fluid_area(optProperties.lookup("fluid_area"));
Switch test_area(optProperties.lookup("test_area"));

dimensionedVector omega_vector(optProperties.lookup("omega_vector")); //转速，rpm

Info<< "Reading pseudo-density field\n" << endl; //创建伪密度场xp——PDE过滤后的
volScalarField xp
(
    IOobject
    (
        "xp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),                                                                  // 伪密度场初始化为voluse
    zeroGradientFvPatchScalarField::typeName      // 密度场边界为零梯度
);
volScalarField xh                                                                   //创建伪密度场xh——Heviside过滤后的
(
    IOobject
    (
        "xh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),                                                                  // 伪密度场初始化为voluse
    zeroGradientFvPatchScalarField::typeName     // 密度场边界为零梯度
);
volScalarField x(xp);                                                           //创建伪密度场x，初始化为xp
volScalarField drho(x);                                                      //创建标量场drho，用于“Heaviside过滤的链式求导”


/* 定义固体区域、流体区域和用于计算直管能量耗散的test_area */
labelList cells_solid,cells_fluid,cells_test;               // 创建列表，用于存储域内的单元序号
if(solid_area)                                                                         // solid_area是个bool型变量，yes/no 
{
  word zoneName="solid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_solid=mesh.cellZones()[zoneID];
  //相当于labelList& cells _solid= mesh.cellZones()[mesh.cellZones().findZoneID("zoneName")];
  forAll(cells_solid, i)
  {
     x[cells_solid[i]]=0;                                                           // 给所有solid单元伪密度x赋值为0
  }
}
if(fluid_area)
{
  word zoneName="fluid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_fluid=mesh.cellZones()[zoneID];
  forAll(cells_fluid, i)
  {
     x[cells_fluid[i]]=1.0;                                                          // 给所有fluid单元伪密度x赋值为1
  }
}
if(test_area)
{
  word zoneName="test_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_test=mesh.cellZones()[zoneID];
  forAll(cells_test, i)
  {
     x[cells_test[i]]=1.0;                                                         // 给所有test单元伪密度x和xh都赋值为1
     xh[cells_test[i]]=1.0;
  }
}
volScalarField gsensVol                                                    // 定义体标量场——体积约束对x的灵敏度 ，初始化为xh
(
    IOobject
    (
        "gsensVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    xh,
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsenshVol(xh);                                          // 定义体标量场——体积约束对xh的灵敏度 

/* 查找边界处的网格 */ 
// find the cells on the boundary
dictionary conFunc = mesh.solutionDict().subDict("constraintFunctionDict_flow");
int nObjPatch =conFunc.lookupOrDefault<scalar>("numberConstraintPatches",2);
wordList conPatchNames = conFunc.lookup("constraintPatchesNames");
label conPatchList[nObjPatch];//patch list，conPatchList在“costunction”计算能量耗散用到了
int iLoop;
for(iLoop=0; iLoop<nObjPatch; iLoop++)
{
    conPatchList[iLoop]=mesh.boundaryMesh().findPatchID(conPatchNames[iLoop]);
}   
