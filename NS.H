for(i=0;i<40;i++)                                                               //相当于while（simple.loop()）,执行40次simple循环
{
MRF.correctBoundaryVelocity(U);
tmp<fvVectorMatrix> tUEqn
(
    fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
  //- fvm::laplacian(nu,U)
  + fvm::Sp(alpha, U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);  //对速度方程进行限制，保证fvOptions功能的顺利实现

if (simple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);//修正速度值，保证fvOptions功能的顺利实现
}
//****************************************//
volScalarField rAU(1.0/UEqn.A());   //A矩阵的逆矩阵
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));//计算A的逆乘以H，也就是H除A
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
MRF.makeRelative(phiHbyA);// 使绝对通量转化为相对通量，因为MRF要用到绝对速度向相对速度的转换。
adjustPhi(phiHbyA, U, p);//调整通量使其满足连续性

tmp<volScalarField> rAtU(rAU);

if (simple.consistent())
{
    rAtU = 1.0/(1.0/rAU - UEqn.H1());
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

tUEqn.clear();//由于simple算法迭代中系数矩阵不断更新，此书将tueqn清空，为修正速度做准备

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (simple.correctNonOrthogonal())//非正交很严重时才用，且不超过2
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (simple.finalNonOrthogonalIter())计算到非正交修正的最后一部时按照3-188计算通量，
    {
        phi = phiHbyA - pEqn.flux();          // fvc::flux(U)  // 通量的值， 该句等同于fvc::interpolate(U) & mesh.Sf()
    }
}

#include "continuityErrs.H"  //计算连续性误差并在屏幕上输出

// Explicitly relax pressure for momentum corrector
p.relax();//显示修正压力

// Momentum corrector
U = HbyA - rAtU()*fvc::grad(p);//采用新的压力修正速度
U.correctBoundaryConditions();//修正边界速度
fvOptions.correct(U);//当使用fvOption时进一步修正速度

U.storePrevIter();
p.storePrevIter();
phi.storePrevIter();

}



