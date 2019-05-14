# Matlab-Picoscope
Script to aquire from a *Picoscope 4244* with *Matlab*.
The script is based on _PS4000_IC_Generic_Driver_Streaming_.

## Instruction ##
Download and set all the Picoscope SDK. Make sure you link the SDK folder to your Matlab path.
Add the script to the SDK folder and work from there.
Connect your Picoscope and set all the parameters in the script (time interval, channels ranges, ...). 
Launch the script and wait that it reaches the point to start the aquisition. When you want launch the aquistion.
Use _CollectData_ to convert the data into vectors and then _int2float_ to convert these vectors to float values for the next analysis.

## Note ##
Most of the modifications from the original script are marked with *Alberto* o *Alb*.

# Authors #
Alberto Battistel and Fabio La Mantia
