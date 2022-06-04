Info<< "Reading TransportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
/* 读取流体属性 */
dimensionedScalar nu(transportProperties.lookup("nu"));                                     //流体粘度
dimensionedScalar alphaMax(transportProperties.lookup("alphaMax"));        //初始流阻系数的最大值 （刚开始优化时放松）
dimensionedScalar alphaMAX(transportProperties.lookup("alphaMAX"));      // 流阻系数的最大值alpha_max
dimensionedScalar alphamax(transportProperties.lookup("alphaMax"));
scalar qu(0.01);                                                                                                                               // 定义插值常数q为0.01
volScalarField alpha(alphaMax*qu*(1-xh)/(qu+xh));                                                      // 体标量场流阻系数插值alpha_f

#include "createMRF.H"                                                      // simpleFoam有
Info<< "Reading field p\n" << endl;
volScalarField p                                                                       // 定义压力场p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field U\n" << endl;                                 // 定义速度场U
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
<<<<<<< HEAD
#include "createPhi.H"

/*设置压力参考*/
=======
#include "createPhi.H"                                                          // 定义通量场phi，面心标量场，48行到57是simpleFoam有的
>>>>>>> 68d448599c810c387334e453522b53f37b270c14
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);//表示从/system/fvsolution/中的SIMPLE部分设置的pRefCell, pRefValue，没有设置采用默认值。

mesh.setFluxRequired(p.name());
<<<<<<< HEAD


singlePhaseTransportModel laminarTransport(U, phi);//单相输运模型中的粘度模型，此处仅为构造粘度模型，为了让用户在算例的/constant/transportantProperties中设置transportModel
/*autoPtr是openFOAM的只能指针，格式为autoptr<类名称>对象名 */
autoPtr<incompressible::turbulenceModel> turbulence
=======
singlePhaseTransportModel laminarTransport(U, phi);
autoPtr<incompressible::turbulenceModel> turbulence// turbulence是个autoPtr
>>>>>>> 68d448599c810c387334e453522b53f37b270c14
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)//表示指向选定的RANS模型
);
//*************************************************// adj-fluid-a
Info<< "Reading field pa\n" << endl;
volScalarField pa                                                                       // 定义伴随压力场pa
(
    IOobject
    (
        "pa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field Ua\n" << endl;
volVectorField Ua                                                                       // 定义伴随速度场Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading/calculating face flux field phi_a\n" << endl;
surfaceScalarField phia                                                              // 定义伴随通量场phi_a，面心标量场
(
    IOobject
    (
        "phi_a",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Ua) & mesh.Sf()                                       // 值为面上的速度u&面矢量sf
);
label paRefCell = 0;                                                                        // pa的参考压力
scalar paRefValue = 0.0;
setRefCell
(
    pa,
    simple.dict(),
    paRefCell,
    paRefValue
);
mesh.setFluxRequired(pa.name());
//****************************************************************//adj-fluid-b
Info<< "Reading field pb\n" << endl;
volScalarField pb
(
    IOobject
    (
        "pb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field Ub\n" << endl;
volVectorField Ub
(
    IOobject
    (
        "Ub",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading/calculating face flux field phi_b\n" << endl;
surfaceScalarField phib
(
    IOobject
    (
        "phi_b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Ub) & mesh.Sf()
);
label pbRefCell = 0;
scalar pbRefValue = 0.0;
setRefCell
(
    pb,
    simple.dict(),
    pbRefCell,
    pbRefValue
);
mesh.setFluxRequired(pb.name());

//*****************************************************//
volScalarField gsensPowerDiss                                                                               // 定义能量耗散约束对x的灵敏度，
(
    IOobject
    (
        "gsensPowerDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ua),                                    // 应为R^u(动量方程)中只有alpha_f*U项与xh有关，所以dJ/dxh=
                                                                                                                                             // ua&d(alpha_f*U)/dxh
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsenshPowerDiss(gsensPowerDiss);                                   // 定义能量耗散约束对xh的灵敏度，
