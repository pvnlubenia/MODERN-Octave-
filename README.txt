==========================================================================================

  MODERN: Mass action Ordinary Differential Equation system to chemical Reaction Network

==========================================================================================

GNU Octave (https://www.gnu.org/software/octave/index) was used to develop the function used here.



========
Function
========The function ode_to_crn returns the chemical reaction network corresponding to a system of ordinary differential equations (ODEs) with mass action kinetics. If the system does not satisfy the Hars-Toth criterion, a message appears stating what needs to be added in a flux.

The output variable 'ode' allows the user to view the complete system of ODEs with all the species and fluxes listed in the 'species' and 'flux' fields, respectively.



=================================
How to fill out 'ode' structure
=================================

'ode' is the input for the function ode_to_crn. It is a structure, representing the system of ODEs, with the following fields:

   - id: name of the model
   - species: a list of all species in the system; this is left blank since incorporated into the function is a step which compiles all species used in the ODEs
   - flux: a list of all fluxes in the system; this is left blank since incorporated into the function is a step which compiles all fluxes used in the ODEs
   - kinetics: a list of all kinetics in the system, each with the following subfields:
        - id: has the following further subfields:
             - var: a string representing the flux number of the kinetics
             - eq: a string representing the kinetic equation
        - species: a list of strings representing the species in the flux
        - order: a list of numbers representing the kinetic order of each species in the flux in the left to right direction (listed  in the same order of the species)
   - equation: a list of all ODEs in the system, each with the following subfields:
        - id: has the following further subfields:
             - var: a string representing the dependent variable (state variable) of the ODE
             - eq: a string representing the ODE flux balance equation
        - flux: a list of fluxes representing the fluxes in the ODE
        - coeff: a list of numbers representing the coefficient of each flux in the ODE in the left to right direction (listed  in the same order of the fluxes)

Note that it is assumed that the ODE system has mass action kinetics.



========
Examples
========

4 examples are included in this folder, all of which came from [1].



===================
Contact Information
===================

For questions, comments, and suggestions, feel free to contact me at pvnlubenia@yahoo.co.uk.


- Patrick Lubenia (19 June 2022)



=========
Reference
=========

   [1] Chellaboina V, Bhat S, Haddad W, Bernstein D (2009) Modeling and analysis of mass-action kinetics. IEEE Control Syst 29(4):60-78. https://doi.org/10.1109/MCS.2009.932926