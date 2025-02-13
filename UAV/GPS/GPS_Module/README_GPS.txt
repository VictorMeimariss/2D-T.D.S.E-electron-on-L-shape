GPS_Module--> Prototype for session launch, we'll be changing settings
and saving data to compare with.

1) Will be setting main GPS Module settings once, then only running with Least Angle Regression (LARS) and model and averaging based on collected data.
2) Setting actual coordinates (mapping on path on google maps), as well as pure GPS coordinates ( walking and riding motorcycle on that path)
to train (LARS) model as well as use averaging to fix randomness errors at each point and then testing.
3) Aligning GPS with API to get access to data at a point and continously in real time.


4) !!! Update GPS coordinates thankfully are fine, maybe no lars needed.