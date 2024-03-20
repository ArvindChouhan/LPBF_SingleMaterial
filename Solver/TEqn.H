int loopcount = 0;
scalar conv = 1;
label LfMaxIter = 10;
scalar LfConverge = 0.001;
alphaT = alpha1;

volScalarField lfPRELOOP = lf;
radiation->correct();

volScalarField factor2 = (scalar(2.0)*rhoCp)/(rho1*(cps*(scalar(1)-lf)+cpl*lf)+rho2*cp2);  
volScalarField loss1 = h1 * (T - Tref) + E*sig*(pow(T,scalar(4.0)) - pow(Tref,scalar(4.0))) ;
volScalarField loss2 = 0.82 * ((L_v * M) / (sqrt(2 * 3.14 * M * R * T))) * P_0 * exp(( L_v * M * (T - T_v))/( R * T * T_v));

while (conv > LfConverge) 
{
	loopcount++;
	volScalarField lfold = lf;
	dimensionedScalar dFdT = 1 / (Tl - Ts);
	volScalarField invF = lfold * (Tl - Ts) + Ts;
	volScalarField dH = alpha1*rho * L * pos(Tl - T) * pos(T - Ts); //alpha1 multiplication
	Sp = - dH/runTime.deltaT() * dFdT;
	Sc = dH/runTime.deltaT() * (lfPRELOOP - lf + dFdT * invF);

	fvScalarMatrix TEqn
        (
             fvm::ddt(rhoCp, T)
          +  fvm::div(rhoCpphi, T)
	  -  fvm::laplacian(k, T)
	  -  fvm::Sp(Sp, T)
          == Sc + radiation->ST(T) - (loss1 + loss2) * mag(fvc::grad(alpha1)) *factor2
       	  //== Sc + (Int - (loss1 + loss2)) * mag(fvc::grad(alpha1)) *factor2
        );

	TEqn.relax();
        fvOptions.constrain(TEqn);
	TEqn.solve();
        fvOptions.correct(T);

	conv = -1;
	lf = (lfold + (dFdT*(T-invF)));
	lf = max(lf,scalar(0));
	lf = min(scalar(1),lf);
	conv = max(mag((lf-lfold)/(lfold+scalar(1e-3)))).value();
	if (loopcount > LfMaxIter) break;
	Info << "\tloopcount: " << loopcount << ", \tRESIDUAL: " << conv << endl;

	loss1 = h1 * (T - Tref) + E*sig*(pow(T,scalar(4.0)) - pow(Tref,scalar(4.0))) ;
	loss2 = 0.82 * ((L_v * M) / (sqrt(2 * 3.14 * M * R * T))) * P_0 * exp(( L_v * M * (T - T_v))/( R * T * T_v));
}

Info << "min/max(T) = " << min(T).value() << ", " << max(T).value() << endl;	

 #include "G_V_CR.H"





