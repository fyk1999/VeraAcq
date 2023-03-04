# fyk_L11_4vDopplerFlashAngles

- Real-time doppler acquistion and display，each doppler image is computed from ne pages of IQ-frame, and each IQ-frame is compounded by na angles of plane wave.

- So far it did not add the Singular-Value-Decomposition(SVD) clutter filtering method, considering the huge computation cost of SVD, it still need a lot of work to realize it in a  real-time ultrasonic machine

# fyk_L11_4vAcquireDopplerIQ_AnglesSC

- Acquire post-compounded IQ data. After one frame in the RcvBuffer is fullfilled, they are reconstructed to pixel-based IQ data and sent to the InterBuffer, and then we use an external function to save this frame of IQ data.

# fyk_L11_4vAcquireDopplerIQ_AnglesSC_AVRF
- The difference of this script with the previous one is that we add some GUI to this. Specifically, the regin where the doopler sequence will acquire is ploted in the B-mode image, and User can adjust this region when the B-mode is displayed.

- A second improvement of this script is that it adopts a method called "Averaging RF", which will transmit two waves in each angle. This method can increase the SNR while it has no need of conducting double DAS opertate.

# fyk_C5_2vAcquireDopplerIQ_AnglesSC
- Adopt the previous script from a linear array to a curved array, so that we can acquire some data of human kidney.