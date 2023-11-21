# BLDCTempMapper
Simulation tool for hobbyists, researchers, etc that provides dynamic temperature estimates. 

Example below uses default values and has several periods of no movement where the current = 0, such as [300s - 400s]. The first third simulates a walking gait of a quadruped robot (based on loose current estimates), second third is a trot, and the last third is a running gait. It can be seen that the rotor and stator have approximately the same temperature with the rotor temp always lagging the stator. In the last image it can be seen that with the defaul motor configured and the current profile, the motor can operature under its max operating temperature in walking and trot, but not during a running gait. 


<img width="834" alt="image" src="https://github.com/dkanh6/BLDCTempMapper/assets/147087124/498aa0bb-c5d4-4fe0-9c6b-952c1bd302bd">

<img width="611" alt="image" src="https://github.com/dkanh6/BLDCTempMapper/assets/147087124/2e094622-ae9d-4d6d-9c09-4b4497f14121">


<img width="689" alt="image" src="https://github.com/dkanh6/BLDCTempMapper/assets/147087124/090d3a11-2b50-48cc-9559-9e826664334b">

<img width="632" alt="image" src="https://github.com/dkanh6/BLDCTempMapper/assets/147087124/2586ab83-eb39-4189-820a-202e580b457b">

