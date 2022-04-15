Info<< "\nslove linear elasticity equation\n" << endl;

initialResidual = 0;
iCorr=0;
do
{
  {
      fvVectorMatrix DEqn
      (
          fvm::d2dt2(D)
       ==
          fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
          	+ ((2.0*M_PI/60.0*omega_vector)^((2.0*M_PI/60.0*omega_vector)^mesh.C())) //omega_vector[1/s] in [r/m]
        + divSigmaExp
        //- fvc::grad(threeKalpha*T)
      );

      initialResidual = DEqn.solve().max().initialResidual();
      if(runTime.outputTime())// 我加的，输出D
        {
          D.write();
        }

      if (!compactNormalStress)
      {
          divSigmaExp = fvc::div(DEqn.flux());
      }
  }

  {
      volTensorField gradD(fvc::grad(D));
      sigmaD = mu*twoSymm(gradD) + (lambda*I)*tr(gradD);

      if (compactNormalStress)
      {
          divSigmaExp = fvc::div
          (
              sigmaD - (2*mu + lambda)*gradD,
              "div(sigmaD)"
          );
      }
      else
      {
          divSigmaExp += fvc::div(sigmaD);
      }
  }

} while (initialResidual > convergenceTolerance && ++iCorr < nCorr);
