void GetRecenteringCorrection(TH2D *hQ, float &xmean, float &ymean)
{
    xmean = hQ->GetMean(1);
    ymean = hQ->GetMean(2);
}

void GetWidthCorrection(TH2D *hQ, float &xdev, float &ydev)
{
    xdev = hQ->GetStdDev(1);
    ydev = hQ->GetStdDev(2);
}

void GetTwistAndRescaleCorrection(TH2D *hQ, float &aplus, float &aminus, float &lambdaplus, float &lambdaminus)
{
    float rho = hQ->GetCorrelationFactor();
    float sigmax = hQ->GetStdDev(1);
    float sigmay = hQ->GetStdDev(2);

    float b = rho*sigmax*sigmay*TMath::Sqrt(2.0*(sigmax*sigmax + sigmay*sigmay
                - 2.0*sigmax*sigmay*TMath::Sqrt(1.0 - rho*rho))/((sigmax*sigmax
                + sigmay*sigmay)*(sigmax*sigmax + sigmay*sigmay)
                + 4.0*(sigmax*sigmay*rho)*(sigmax*sigmay*rho)));

    aplus = TMath::Sqrt(2.0*sigmax*sigmax - b*b);
    aminus = TMath::Sqrt(2.0*sigmay*sigmay - b*b);

    lambdaplus = b/aplus;
    lambdaminus = b/aminus;
}

void DoCorrections(float &qx, float &qy, float *corrections)
{
    // 1) Recentering
    qx -= corrections[0];
    qy -= corrections[1];

    // 2) Twist
    qx = (qx - corrections[7]*qy)/(1.0 - corrections[7]*corrections[7]);
    qy = (qy - corrections[6]*qx)/(1.0 - corrections[6]*corrections[6]);

    // 3) Rescale
    if (corrections[4]==0 || corrections[5]==0) {
        std::cout << "Correction warning: aplus or aminus is equal to zero, rescale not applied!" << std::endl;
    } else {
        qx /= corrections[4];
        qy /= corrections[5];
    }
}

float *CalculateCorrections(TH2D *hQ)
{
    static float corrections[8] = {0};
    float xmean, ymean, xdev, ydev, aplus, aminus, lambdaplus, lambdaminus;
    GetRecenteringCorrection(hQ, xmean, ymean);
    GetWidthCorrection(hQ, xdev, ydev);
    GetTwistAndRescaleCorrection(hQ, aplus, aminus, lambdaplus, lambdaminus);
    corrections[0] = xmean;
    corrections[1] = ymean;
    corrections[2] = xdev;
    corrections[3] = ydev;
    corrections[4] = aplus;
    corrections[5] = aminus;
    corrections[6] = lambdaplus;
    corrections[7] = lambdaminus;
    return corrections;
}
