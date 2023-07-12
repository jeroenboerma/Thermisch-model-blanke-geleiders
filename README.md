# Thermisch-model-blanke-geleiders

The model for the determination of the load capacity of bare conductors exists of 6 scripts. 

The first script is called **Properties.py** and defines all parameters and constants needed to define the physical system. Properties.py contains 7 functions with the following functionalities:
1. Component():
   #Define component
    def Component():
      Material = "Copper"     #Choose "Copper" or "Aluminium"
      Shape = "Rectangular"   #Choose "Cylindrical", "Rectangular" or "Stranded"
      Core = "Solid"          #Choose "Solid", "Hollow" or "Steel"
      Location = "Inside"    #Choose "Inside" or "Outside"
      Painted = "No"          #Choose "Yes" of "No"
    return Material, Shape, Core, Location, Painted

Component() is the menu to choose a specific type of conductor. You choose the material, shape, core, location and whether the component is painted or not. Regarding the core, only the shape "Cylindrical" can be hollow, the others are "Solid". It is advised to leave Painted on "No".

2. select_material(): select_material() defines material properties of copper, aluminium and aluminium alloys.
3. component_paramaters(): this function defines the dimensions of the conductor. The length L of the conductor can be chosen but as long as the conductor is not shorter than 1 meter, it does not really matter. This can be left untouched. For all shapes only the following dimensions have to be chosen:
    - Rectangular: Height and Width
    - Cylindrical solid: Outer diameter
    - Cylindrical hollow: Outer diameter and wall thickness
    - Stranded: Number of strands, wire diameter and outer conductor diameter (don't worry, ill provide documentation to find this)
IMPORTANT: Fill in all dimension in [mm], but know that the model works with SI units. For this reason you often see the term: *(1/1000)

4. ambient_parameters(): defines air constants. Only parameter is air velocity. The DIN prescribes v = 0.6 m/s for outside installations, this is a standardized value.
5. PhyC0(): Physics constants
6. Wheather(): defines solar irradiatiance and solar absorption coefficient. Solar irradiance has a big influence on the outcome of the model, 900 W/m2 is average, 1200 W/m2 is maximum for the Netherlands.
7. Rac(): calculates the AC resistance of all shape categories. This should NOT be changed. Mathematical framework can be found in my internship report.

**The second script is Heat_flows.py**
In heat_flows.py, physical formulas are used to define the four dominant heat flows: Ohmic heat generation, Solar irradiance, Radiation and Convection. Each heat flow has its own function. There is a fifth function calling heat capacity properties. Temperatures are always in Kelvin 

1. Qemis(): defines the heat lost by grey body radiation. This is a standard physical formula.
2. Qirr(): If conductor is inside, the irradiance is zero. If conductor is outside, irradiance is defined.
3. Qgen(): Defines Ohmic heat. First the AC resistance is defined by calling the Rac() function from Properties.py.
4. Qconv(): Defines the heat lost by convection. There are two forms of convection: natural and forced. IF Location == "Inside" --> natural convection. IF Location == "Outside" --> forced convection. If preferred to have forced convection in an inside situation(because of ventilation for example), the Location should be set to "Outside" and the solar irradiance should be set to zero in Properties.Wheather()  
For every conductor shape the natural and forced convection are defined seperately. Mathematical framework can be found in the internship report. In this function you will find IF-statements like:
            if Ra > 10**12: #Limit of NuD equation
                print('Rayleigh number bigger than 10^12, outcome not valid!') 
Officially, the physical formulas do not apply when these type of statements are true. It is expected that these statements are never triggered, but if so, check the magnitude of all numbers you have filled in. This may be a reason for exceptional values triggering the IF-statement.
5. Qabsorb(): just calls material properties.

**The third script is Temperature_profile.py**
This script contains two functions.
1. Times(): defines the timeframe for how long you want to run the model. You can set the amount of hours for running a single temperature curve. For most conductors 4 hours is more than enough to let the curve flatten out, but for bigger or solid conductors you might need 6 hours. The timesteps are set at dt = 1 [s]. This time resolution if fine.

2. Temp_profile(I, Te): this is the main function for calculating the heat balance and updating the temperature. First the Times() function is called to know how many iterations the model should run. Then empty arrays are created for all parameters that need to be calculated.
The iterations start at the FOR-loop. In the first iteration the temperature of the conductor T is set to Te, the ambient temperature. This means the conductor will always start at ambient temperature. If necessary, a different value can be filled in, but don't forget to set it back to "Te" for convenience. After setting the initial temperature, all heat flow functions are called and the heat balance is calculated. Then the temperature step is calculated and added to the temperature at the start of the iteration. The FOR-loop starts again with the updated temperature and skips the first two lines as "i" is non-zero.

**The fourth script is Tfinal_and_tau.py**
In this script a function is defined which determines the final temperature and the time constant of the modelled warm-up curve. The first FOR-loop finds the iteration at which the temperature is within 0.01 of the final temperature. This is done because this iteration number is needed to define the time constant correctly. One odd thing is that the function which calculates the complete temperature curve is only called in this script in line 20. It would have been more logical to do this in a main script, so this would be an improvement of the model. If ever searching for the calling of the temperature curve, then this is where you should look.

**The fifth script is I_vs_Te.py**
At the top of this script, after calling libraries and functions, the ambient temperature "Te" and current "I" are called. This script plots the warm-up curve for a constant temperature and constant current. There are THREE possibilities:
- select single ambient temperature and single current. Result is a single curve.
- select single ambient temperture and array of currents. Result is a plot with a curve per current.
- select array of ambient temperature and single current. Result is a plot with a curve per ambient temperature.

After selecting the ambient temperature and current empty arrays are created to fill in the for loop. The for-loop iterates over the array of Te or I and calls the warm-up curve for every combination of Te and I. The warm-up curves are plotted, and if there is more than one current value assigned, the dependency between the current and final temperature is plotted. The last plot shows the heat flows. "k" refers to the curve number if multple curves are plotted.

IMPORTANT: This script can mainly be used to check whether the outcome of the model corresponds with values listed in load capacity tables. To determine the new load capacity it is better to use the next script.

**The sixth script is Det_Imax.py**
In this script you have THREE possibilities:
- select single ambient temperature and single current. Result is a single curve.
- select single ambient temperture and array of currents. Result is a plot with a curve per current.
- select array of ambient temperature and single current. Result is a plot with a curve per ambient temperature.

Also, fill in the nominal current as listed by the ABB pocket book. 
This script outputs load capacity values expressed either in ampere or percentage of the nominal current. To choose which one you want, put "Yes" or "No" behind I_want_the_load_capacity_in_percent_of_Inom.
The script calculates the load capacity based on an heat balance, this happens inside the for-loop. After the loop, the relevant temperatures are set from Kelvin to degrees Celsius.
Depending on which combination of ambient temperature and current the model outputs either a graph showing the relation between ambient temperature or maximum temperature versus the load capacity, or prints out a single value.

**General remark**
Script 1,2,3 and 4 are essential to run either script 5 or 6. 5 and 6 themselves are independent of each other. My advise is to always have all 6 scripts open in consecutive order. When using the model, you should only need to change parameters in script 1 Properties.py and in script 5 or 6 I_vs_Te.py and Det_Imax.py. Only Temperature_profile.py is used to change the run time. 
