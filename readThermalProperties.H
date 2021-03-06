Info<< "Reading thermal properties\n" << endl;
IOdictionary thermalProperties
(
    IOobject
    (
        "thermalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar rhoc(thermalProperties.lookup("rhoc"));                   // 流体密度乘以比热容
dimensionedScalar kf(thermalProperties.lookup("kf"));                               // 流体热传导系数
dimensionedScalar ks(thermalProperties.lookup("ks"));                               // 固体热传导系数
dimensionedScalar Q(thermalProperties.lookup("Q"));                                 // 体热源 Q/(rho*Cp)
//dimensionedScalar Q("Q", dimensionSet(0,0,-1,1,0,0,0),100000000.0/4180000);

volScalarField DT                                                                                                               // 体标量场热传导系数
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    (ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoc                                                                     // 热传导系数k插值
);
Info<< "Reading field T\n" << endl;
    volScalarField T                                                                       // 定义温度场T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
Info<< "Reading field Tb\n" << endl;
   volScalarField Tb                                                                       // 定义伴随温度场Tb
   (
       IOobject
       (
           "Tb",
           runTime.timeName(),
           mesh,
           IOobject::MUST_READ,
           IOobject::AUTO_WRITE
       ),
       mesh
   ); 
volScalarField T2(T);                                                                    // T2 应该没有用到
volScalarField fsensMeanT                                                       // 定义体标量场fsensMeanT，平均温度对x的灵敏度
(
    IOobject
    (
        "fsensMeanT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ub),   // 平均温度的灵敏度初始化
    zeroGradientFvPatchScalarField::typeName
);
volScalarField fsenshMeanT(fsensMeanT);                      // 定义体标量场fsenshMeanT，平均温度对xh的灵敏度
