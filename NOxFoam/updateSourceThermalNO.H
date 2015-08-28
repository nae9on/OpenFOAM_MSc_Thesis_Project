
// Calculating the source term
volScalarField num1 ("num1", scalar(1) - (kr1*kr2*concNO*concNO)/(kf1*kf2*(concN2+smallConc)*(concO2+smallConc)) );
//volScalarField num1 ("num1", concN2*concO2 - (kr1*kr2)*concNO*concNO/(kf1*kf2) );
//Info<< "num computed\n" << endl;

volScalarField denom1 ("denom1", scalar(1) + (kr1*concNO)/(kf2*(concO2+smallConc)+kf3*concOH) );
//volScalarField denom1 ("denom1", concN2*concO2 + (kr1*concNO*concN2*concO2)/(kf2*concO2+kf3*concOH) );
//Info<< "denom computed\n" << endl;

sourceThermalNO = M_NO * scalar(2) * kf1 * concO * concN2 * num1 / denom1;

//Info<< "source computed\n" << endl;

//Info<< "Source Thermal NO Updated .................................................................................\n" << endl;
