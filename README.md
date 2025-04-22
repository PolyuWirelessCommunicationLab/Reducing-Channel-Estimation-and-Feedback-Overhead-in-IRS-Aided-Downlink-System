# Reducing-Channel-Estimation-and-Feedback-Overhead-in-IRS-Aided-Downlink-System
This code is for paper: [R. Wang, Z. Wang, L. Liu, S. Zhang, and S. Jin, “Reducing channel estimation and feedback overhead in IRS-aided downlink system: A quantize-then-estimate approach,” IEEE Trans. Wireless Commun., vol. 24, no. 2, pp. 1325–1338, Feb. 2025.](https://ieeexplore.ieee.org/abstract/document/10786374)
# Abstract
Channel state information (CSI) acquisition is essential for the base station (BS) to fully reap the beamforming gain in intelligent reflecting surface (IRS)-aided downlink communication systems. Recently, [1] revealed a strong correlation in different users’ cascaded channels stemming from their common BS-IRS channel component, and leveraged such a correlation to significantly reduce the pilot transmission overhead in IRSaided uplink communication. In this paper, we aim to exploit the above channel property to reduce the overhead for both pilot and feedback transmission in IRS-aided downlink communication. Note that in the downlink, the distributed users merely receive the pilot signals containing their own CSI and cannot leverage the correlation in different users’ channels, which is in sharp contrast to the uplink counterpart considered in [1]. To tackle this challenge, this paper proposes a novel “quantize-then-estimate” protocol in frequency division duplex (FDD) IRS-aided downlink communication. Specifically, the users quantize and feed back their received pilot signals, instead of the estimated channels, to the BS. After de-quantizing the pilot signals received by all the users, the BS estimates all the cascaded channels by leveraging their correlation, similar to the uplink scenario. Under this protocol, we manage to propose efficient user-side quantization and BS-side channel estimation methods. Moreover, we analytically quantify the pilot and feedback transmission overhead to reveal the significant performance gain of our proposed scheme over the conventional “estimate-then-quantize” scheme.
# Citation
@article{Wang2025reducing,
  author={Wang, Rui and Wang, Zhaorui and Liu, Liang and Zhang, Shuowen and Jin, Shi},<br> 
  journal={IEEE Trans. Wireless Commun.},<br> 
  title={Reducing Channel Estimation and Feedback Overhead in IRS-Aided Downlink System: A Quantize-Then-Estimate Approach}, <br>
  year={2025},<br>
  volume={24},<br>
  number={2},<br>
  pages={1325-1338},<br>
  keywords={Channel estimation;Downlink;Protocols;Quantization (signal);Correlation;Uplink;Feeds;Vectors;Wireless communication;Transmitting antennas;Intelligent reflecting surface (IRS);channel 
  estimation;channel feedback;distributed source coding},<br>
}
# Note
The code is provided for the benefit of better understanding the results, and is not meant to be used in production.
