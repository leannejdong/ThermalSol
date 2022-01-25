# Thermal Simulations of Networks

As a continuation of the hydraulic simulations of the verona network in [HydraulicNetwork](https://github.com/leannejdong/HydraulicNetwork),
we perform thermal simulations in this project.

Our code first obtain the required mass flow data from the outputs of [HydraulicNetwork](https://github.com/leannejdong/HydraulicNetwork),
then in each step corrects the direction of each pipe (if the mass flow rate is negative we reverse the input and output nodes)
It starts solving the equation for all pipes and if one is connected to different pipes, the final temperature for that particular node become
the average of the output temperatures.