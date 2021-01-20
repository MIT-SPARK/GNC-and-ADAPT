function [weights, th1, th2] = gncWeightsUpdate(weights, mu, residuals, barc2)
%GNC - GNC weights update for the TLS cost function
% Reference:
%    - "Graduated Non-Convexity for Robust Spatial Perception: From Non-Minimal Solvers to Global Outlier Rejection"
%      Yang, Antonante, Tzoumas, Carlone (2020). IEEE Robotics and Automation Letters (RA-L), 5(2), 1127–1134.
%      https://arxiv.org/pdf/1909.08605.pdf
%    - "Outlier-Robust Estimation: Hardness, Minimally-Tuned Algorithms, and Applications" 
%      Antonante, Tzoumas, Yang, Carlone (2020).
%      https://arxiv.org/pdf/2007.15109.pdf
%

% Author: Pasquale Antonante
% email: antonap@mit.edu
% Date: 2021-01-06

th1 = (mu+1)/mu * barc2;
th2 = (mu)/(mu+1) * barc2; % th1 > th2
for k = 1:length(residuals)
    if residuals(k) - th1 >= 0
        weights(k) = 0;
    elseif residuals(k) - th2 <= 0
        weights(k) = 1;
    else
        weights(k) = sqrt( barc2*mu*(mu+1)/residuals(k) ) - mu;
        assert(weights(k)>= 0 && weights(k) <=1, 'weights %g calculation wrong!', weights(k));
    end
end
end
