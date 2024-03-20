

///////////////////////////////////////////////////////////////////////////////
if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    tmp<volScalarField> trSubDeltaT;

    if (LTS)
    {
        trSubDeltaT =
            fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
    }

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

rho = alpha1*rho1 + alpha2*rho2;
rhoCp = alpha1*rho1*(cps*(scalar(1)-lf)+cpl*lf) +  alpha2*rho2*cp2;

k =  alpha1*(ks*(scalar(1)-lf)+kl*lf) +  alpha2*k2;

beta = alpha1*beta1 +  alpha2*beta2;

A = -C * alpha1 *sqr(scalar(1)-lf)/(pow(lf,scalar(3))+ b);

sigm = a_sig + b_sig * (T-Tl);
sigm = max(sig_limit,sigm);     //limiting


