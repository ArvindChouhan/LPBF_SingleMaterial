    volVectorField gradAlpha = fvc::grad(alpha1);
    volVectorField nHat = gradAlpha/(mag(gradAlpha) + mixture.deltaN()); 
    volScalarField factor1 = (scalar(2.0)*rho)/(rho1+rho2);
    kappa = mixture.Arv();

    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U) 
      - fvm::Sp(A,U) 
      - sigm * kappa * nHat  * mag(fvc::grad(alpha1)) * factor1
      - dsigmadt * (fvc::grad(T) - nHat * (nHat & fvc::grad(T))) * mag(fvc::grad(alpha1)) * factor1 
      - 0.54 * P_0 * exp(( L_v * M * (T - T_v))/( R * T * T_v)) * nHat * mag(fvc::grad(alpha1)) * factor1
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
		
		fvc::interpolate(rho* beta * (T - Tref)) * (g & mesh.Sf())
		+
                (
               //     mixture.surfaceTensionForce()*fvc::interpolate(factor1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
