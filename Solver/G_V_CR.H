
        gradT = fvc::grad(T);
        
        volVectorField gradTn = gradT / (mag(gradT)+GC);                  // Normal unit vector to the surface

        V = gradTn & Vlaser;
        
        CR = mag(gradT) * V ;

    
       















