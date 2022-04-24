# IAPWS-IF97-OpenFOAM-v8
***

IAPWS-IF97 Properties of Water and Steam library compatible for OpenFOAM v8.



- Installation introduction
	 - run following command in terminal to compile the lib
	
	   ```bash
	   wmake libso
	   ```
	
- usage
	- add the compiled lib in file controlDict like that:
	
	  ```c++
	  libs
	  (
	    "libfluidThermophysicalModelsNew.so"
	  )
	  ```
	  
	 - set IAPWS-IF97 in constant/thermophysicalProperties simply with:
	
	   ```c++
	   /*--------------------------------*- C++ -*----------------------------------*\
	     =========                 |
	     \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
	      \\    /   O peration     | Website:  https://openfoam.org
	       \\  /    A nd           | Version:  8
	        \\/     M anipulation  |
	   \*---------------------------------------------------------------------------*/
	   FoamFile
	   {
	       version     2.0;
	       format      ascii;
	       class       dictionary;
	       location    "constant";
	       object      thermophysicalProperties;
	   }
	   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	   
	   thermoType  IAPWSThermo;
	   ```

	- run the case as normal:
	
	  ```c++
	  decomposePar
	  mpirun -np 4 buoyantSimpleFoam -parallel
	  ```
	
- result of test case **buoyantCavity_IAWPS**:

   ![Temperature](buoyantCavity_IAWPS\Temperature.png)![p_rgh](buoyantCavity_IAWPS\p_rgh.png)![velocity-Y](buoyantCavity_IAWPS\\velocity-Y.png)
