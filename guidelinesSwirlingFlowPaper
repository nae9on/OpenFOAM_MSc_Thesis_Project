Stabilization of Non-premixed Turbulent Combustion by a Swirling Air Flow

http://www.scialert.net/fulltext/?doi=ajaps.2012.445.459&org=12



In the present case, the fuel is a liquid (fuel oil n°2).

Geometry and Mesh:
//////////////////////// 
In this study, the combustion chamber geometry is a typical case of those used in industrial drying furnaces.

To obtain a non-premixed turbulent flame, the fuel oil n°2 is injected through the internal tube of the burner.
At the same time, primary air at ambient temperature is introduced through the annular space.
Finally, secondary air is introduced through 24 circular orifices on the edges of the front of our firebox.

Mesh: It is important to note that OpenFOAM is strictly a 3D code. To find an initial solution for a 2D axisymmetric problem, the mesh must be 3D first with one cell thick and having no solution in Z-direction. We must ensure that the opening of the thickness of both sides between the Y-axis does not exceed the angle whose apex would be on the Y-axis and between 2.5 and 5°.

Boundary Conditions
////////////////////////
Boundary conditions:  
fuel inlet velocity=2.34 m sec-1, temperature=120 degree celcius
primary air inlet velocity=15 m sec-1, temperature=17 degree celcius
The outlet conditions are related to the pressure. 
Temperature of the nozzle was fixed to 980°C 
secondary air, the inlet velocity is a rotating profile field whose variation depends on the ratio of flow rates injection. 
About the walls, they are supposed to be adiabatic.

//MyNote:
You can define an adiabatic wall by setting a zero heat flux condition. This is the default condition for all walls.

Turbulence
////////////////////////
The model chosen is RANS (Reynolds Averaged Navier-Stokes) k-ε standard because of its robustness, accuracy, low cost and rich documentation. Its weakness in walls zone can be offset by the Standard Model Wall Function (SWF).

Combustion
////////////////////////
For combustion gas, there are three types of solvers all operating in unsteady state(transient models):
reactingFOAM (Lundstrom, 2008), Xoodles and XiFoam. Only reactingFOAM models the non-premixed turbulent combustion. The concept it uses for the chemical species is the PaSR (Partially Stirred Reactor) which is a modified version of the EDC (Eddy Dissipation Concept) where the chemical time scale is handled differently:

NB: the reaction mechanism is imported from Chemkin (software tool for solving complex chemical kinetics problems):

• 	ReactingFOAM Code:
• 	Calculate chemical reaction based on turbulent and chemical timescales
• 	Calculate of the density
• 	Calculate velocity/pressure fields
• 	Read species and feed them to the chemistry solver using Chemkin table
• 	Calculate the temperature from the chemical reactions enthalpy lookup
• 	Calculate the pressure field using PISO
• 	Correct the turbulence (pressure-corrector)
• 	Update the density from the temperature
• 	Return to step 1

Temporal discretisation: To avoid instability due to the simultaneous calculation of turbulence and combustion, the first solver to use is simpleFOAM which is a steady-state incompressible turbulence model. This is done to calculate the cold flow which will be considered as our initial solution. Knowing that reactingFOAM is an unsteady solver, the time step depends strongly on the Courant number Cr. This is possible only with the CFL (Courant Friedrichs Lewy) condition:

If 0<Cr = 0.5, we have more stability and less speed. If 1>Cr = 0.5, we have more speed and less stability. In this case, Cr = 0.2.

Interpolation schemes:
• 	Pressure: LimitedLinear 1 (second order bounded scheme)
• 	Velocity: LimitedLinear V (TVD scheme recommended for swirl)
• 	Turbulence: Upwind (first order bounded scheme)
• 	Species: Upwind
• 	Energy: Upwind
• 	Pressure-velocity coupling: PISO (Pressure Implicit with Splitting of Operators)

Stabilization of combustion: This is done by increasing the residence time of the flame and creating recirculation in the reaction zone. The effect produced is to favor the mixture and give to the flame a more compact form (Poireault, 1997) (Fig. 9-13). The axial evolution of swirling temperature profile is almost linear (Fig. 11). It means that swirl configuration provides better and easiest control of temperature in flow direction.

