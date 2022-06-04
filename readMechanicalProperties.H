Info<< "Reading mechanical properties\n" << endl;
IOdictionary mechanicalProperties                                       // 创建mechanicalProperties字典
(
    IOobject
    (
        "mechanicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
int iCorr=0,nCorr=150;                                                                          // 这两行是solidDisplacement的fvSolution里stressAnalysis子字典里的内容
scalar convergenceTolerance(1e-5),initialResidual(0);           // nCorr就是“nCorrectors”（稳态为1），convergenceTolerance就是“D”，iCorr没有用到

Switch planeStress(mechanicalProperties.lookup("planeStress"));// 是否平面应力问题
dimensionedScalar rhoE("rhoE", dimMass/dimLength/sqr(dimTime),readScalar(mechanicalProperties.lookup("rhoE")));// 定义弹性模量及单位
dimensionedScalar rho("rho",  dimMass/dimVolume,readScalar(mechanicalProperties.lookup("rho")));// 定义密度及单位
scalar Po(readScalar(mechanicalProperties.lookup("Po")));           //定义无量刚数泊松比
dimensionedScalar E(rhoE/rho);                                                                   // Normalising E，弹性模量正则化，solidDisplacement有，这里的E就是Emax
dimensionedScalar Emin(E*1e-9);                                                                 // Emin=Emax*10^-9
double qd=0.1;                                                                                                         // 插值常数q=0.1
volScalarField mu(qd*(1-xh)/(qd+xh)*(E-Emin)/(2.0*(1.0+Po))+Emin/(2.0*(1.0+Po)));//拉梅常数插值,mu=E/(2.0*(1.0 + nu))
volScalarField lambda(qd*(1-xh)/(qd+xh)*Po*(E-Emin)/((1.0+Po)*(1.0-2.0*Po))+Po*Emin/((1.0+Po)*(1.0-2.0*Po)));//lambda=nu*E/((1.0 + nu)*(1.0 - 2.0*nu))

Info<< "Reading field D\n" << endl;
volVectorField D                                                                                                        // 定义位移场
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField gradD                                                                                            // 定义位移梯度张量场
(
    IOobject
    (
        "gradD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
      fvc::grad(D)                                                                                                          // 显示求梯度
);

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD                                                                           // 定义应力场，对陈张量场
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mu*twoSymm(gradD) + (lambda*I)*tr(gradD)                                       // 应力应变公式，grad(D)就是应变，twoSymm()Return twice the symmetric part of a tensor
);

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp                                                                              // divSigmaExp就是 fvc::div(sigmaD)，显示计算sigma的散度
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::div(sigmaD)
);

Switch compactNormalStress("yes");                                                          // 案例也是默认yes

if (compactNormalStress)                                                                                   // 83-91，原有
{
    divSigmaExp -= fvc::laplacian(2*mu + lambda, D, "laplacian(DD,D)");
}
else
{
    divSigmaExp -= fvc::div((2*mu + lambda)*fvc::grad(D), "div(sigmaD)");
}
mesh.setFluxRequired(D.name());

volScalarField gsensCompliance                                                                      // 定义柔度对x的敏度
(
    IOobject
    (
        "gsensCompliance",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -gradD&&(rho*(3*xh*xh*(E-Emin)/(2.0*(1.0 + Po))*twoSymm(gradD))),// 敏度公式
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsenshCompliance(gsensCompliance);                           // 定义柔度对xh的敏度


