(* Created with the Wolfram Language for Students - Personal Use Only : www.wolfram.com *)
{ms2, 6*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
    (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 + 
    (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx202 + 
    (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
    (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
    (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx212 + 
    (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220 + 
    (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 + 
    (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222) + 
  6*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + Kappa02*mDxbar220) + 
    Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + Kappa12*mDxbar220) + 
    Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + Kappa22*mDxbar220) + 
    Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
    Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
    Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
    Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + Kappa02*mDxbar222) + 
    Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + Kappa12*mDxbar222) + 
    Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + Kappa22*mDxbar222)) + 
  4*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
    Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + 
    Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
    Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211)) + 
  4*((Lambda1200^2 + Lambda1201^2)*mH2I200 + 
    (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
    (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + 
    (Lambda1210^2 + Lambda1211^2)*mH2I211) + 
  6*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + 
    Kappa20^2 + Kappa21^2 + Kappa22^2)*ms2 + 
  4*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*ms2 + 
  4*Lambdax^2*(mHd2 + mHu2 + ms2) - 8*gN^2*MassBp^2*QSp^2 + 
  2*gN^2*QSp*(3*(md200 + md211 + md222)*Qdp + 
    3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
    3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
    2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
    2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
    2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + 
    (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
  6*(TKappa00^2 + TKappa01^2 + TKappa02^2 + TKappa10^2 + TKappa11^2 + 
    TKappa12^2 + TKappa20^2 + TKappa21^2 + TKappa22^2) + 
  4*(TLambda1200^2 + TLambda1201^2 + TLambda1210^2 + TLambda1211^2) + 
  4*TLambdax^2, 
 (-4*(-4*g1^2*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
      Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*MassB^2 - 
    6*g1^2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*
     MassB^2 - 80*g3^2*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*MassG^2 - 
    30*g2^2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*
     MassWB^2 - 2*g1^2*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx202 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
      (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx212 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 + 
      (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222) - 
    40*g3^2*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx202 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
      (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx212 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 + 
      (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222) + 
    15*((Kappa00*(Kappa00*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
        Kappa01*(Kappa01*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
        Kappa02*(Kappa02*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)))*
       mDx200 + (Kappa00*(Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa10*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
        Kappa01*(Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa11*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
        Kappa02*(Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa12*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)))*
       mDx201 + (Kappa00*(Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
            Kappa02*Kappa22) + Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
            Kappa12*Kappa22) + Kappa20*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
        Kappa01*(Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
            Kappa02*Kappa22) + Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
            Kappa12*Kappa22) + Kappa21*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
        Kappa02*(Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
            Kappa02*Kappa22) + Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
            Kappa12*Kappa22) + Kappa22*(Kappa20^2 + Kappa21^2 + Kappa22^2)))*
       mDx202 + (Kappa10*(Kappa00*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
        Kappa11*(Kappa01*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
        Kappa12*(Kappa02*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)))*
       mDx210 + (Kappa10*(Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa10*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
        Kappa11*(Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa11*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
        Kappa12*(Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa12*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)))*
       mDx211 + (Kappa10*(Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
            Kappa02*Kappa22) + Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
            Kappa12*Kappa22) + Kappa20*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
        Kappa11*(Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
            Kappa02*Kappa22) + Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
            Kappa12*Kappa22) + Kappa21*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
        Kappa12*(Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
            Kappa02*Kappa22) + Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
            Kappa12*Kappa22) + Kappa22*(Kappa20^2 + Kappa21^2 + Kappa22^2)))*
       mDx212 + (Kappa20*(Kappa00*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
        Kappa21*(Kappa01*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
        Kappa22*(Kappa02*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)))*
       mDx220 + (Kappa20*(Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa10*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
        Kappa21*(Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa11*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
        Kappa22*(Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa12*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)))*
       mDx221 + (Kappa20*(Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
            Kappa02*Kappa22) + Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
            Kappa12*Kappa22) + Kappa20*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
        Kappa21*(Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
            Kappa02*Kappa22) + Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
            Kappa12*Kappa22) + Kappa21*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
        Kappa22*(Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
            Kappa02*Kappa22) + Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
            Kappa12*Kappa22) + Kappa22*(Kappa20^2 + Kappa21^2 + Kappa22^2)))*
       mDx222) + 
    15*(Kappa00*(Kappa00*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
          (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
          (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220) + 
        Kappa10*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx201 + 
          (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx211 + 
          (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx221) + 
        Kappa20*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx202 + 
          (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx212 + 
          (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx222)) + 
      Kappa01*(Kappa01*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
          (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
          (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220) + 
        Kappa11*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx201 + 
          (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx211 + 
          (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx221) + 
        Kappa21*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx202 + 
          (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx212 + 
          (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx222)) + 
      Kappa02*(Kappa02*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
          (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
          (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220) + 
        Kappa12*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx201 + 
          (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx211 + 
          (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx221) + 
        Kappa22*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx202 + 
          (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx212 + 
          (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx222)) + 
      Kappa10*(Kappa00*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
           mDx200 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx210 + 
          (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx220) + 
        Kappa10*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
           mDx201 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
          (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221) + 
        Kappa20*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
           mDx202 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx212 + 
          (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx222)) + 
      Kappa11*(Kappa01*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
           mDx200 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx210 + 
          (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx220) + 
        Kappa11*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
           mDx201 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
          (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221) + 
        Kappa21*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
           mDx202 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx212 + 
          (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx222)) + 
      Kappa12*(Kappa02*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
           mDx200 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx210 + 
          (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx220) + 
        Kappa12*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
           mDx201 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
          (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221) + 
        Kappa22*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
           mDx202 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx212 + 
          (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx222)) + 
      Kappa20*(Kappa00*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
           mDx200 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
           mDx210 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx220) + 
        Kappa10*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
           mDx201 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
           mDx211 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx221) + 
        Kappa20*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
           mDx202 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
           mDx212 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222)) + 
      Kappa21*(Kappa01*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
           mDx200 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
           mDx210 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx220) + 
        Kappa11*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
           mDx201 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
           mDx211 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx221) + 
        Kappa21*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
           mDx202 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
           mDx212 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222)) + 
      Kappa22*(Kappa02*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
           mDx200 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
           mDx210 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx220) + 
        Kappa12*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
           mDx201 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
           mDx211 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx221) + 
        Kappa22*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
           mDx202 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
           mDx212 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222))) - 
    2*g1^2*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
        Kappa02*mDxbar220) + Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
        Kappa12*mDxbar220) + Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
        Kappa22*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + 
        Kappa02*mDxbar221) + Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + 
        Kappa12*mDxbar221) + Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + 
        Kappa22*mDxbar221) + Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
        Kappa02*mDxbar222) + Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
        Kappa12*mDxbar222) + Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
        Kappa22*mDxbar222)) - 
    40*g3^2*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
        Kappa02*mDxbar220) + Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
        Kappa12*mDxbar220) + Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
        Kappa22*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + 
        Kappa02*mDxbar221) + Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + 
        Kappa12*mDxbar221) + Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + 
        Kappa22*mDxbar221) + Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
        Kappa02*mDxbar222) + Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
        Kappa12*mDxbar222) + Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
        Kappa22*mDxbar222)) + 
    15*(Kappa00*((Kappa00*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22))*
         mDxbar200 + (Kappa01*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22))*
         mDxbar210 + (Kappa02*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22))*
         mDxbar220) + Kappa10*
       ((Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa10*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22))*
         mDxbar200 + (Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa11*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22))*
         mDxbar210 + (Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa12*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22))*
         mDxbar220) + Kappa20*
       ((Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
          Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
          Kappa20*(Kappa20^2 + Kappa21^2 + Kappa22^2))*mDxbar200 + 
        (Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
          Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
          Kappa21*(Kappa20^2 + Kappa21^2 + Kappa22^2))*mDxbar210 + 
        (Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
          Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
          Kappa22*(Kappa20^2 + Kappa21^2 + Kappa22^2))*mDxbar220) + 
      Kappa01*((Kappa00*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22))*
         mDxbar201 + (Kappa01*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22))*
         mDxbar211 + (Kappa02*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22))*
         mDxbar221) + Kappa11*
       ((Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa10*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22))*
         mDxbar201 + (Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa11*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22))*
         mDxbar211 + (Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa12*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22))*
         mDxbar221) + Kappa21*
       ((Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
          Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
          Kappa20*(Kappa20^2 + Kappa21^2 + Kappa22^2))*mDxbar201 + 
        (Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
          Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
          Kappa21*(Kappa20^2 + Kappa21^2 + Kappa22^2))*mDxbar211 + 
        (Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
          Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
          Kappa22*(Kappa20^2 + Kappa21^2 + Kappa22^2))*mDxbar221) + 
      Kappa02*((Kappa00*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22))*
         mDxbar202 + (Kappa01*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22))*
         mDxbar212 + (Kappa02*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
          Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22))*
         mDxbar222) + Kappa12*
       ((Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
          Kappa10*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22))*
         mDxbar202 + (Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa11*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22))*
         mDxbar212 + (Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
            Kappa02*Kappa12) + Kappa12*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
          Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22))*
         mDxbar222) + Kappa22*
       ((Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
          Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
          Kappa20*(Kappa20^2 + Kappa21^2 + Kappa22^2))*mDxbar202 + 
        (Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
          Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
          Kappa21*(Kappa20^2 + Kappa21^2 + Kappa22^2))*mDxbar212 + 
        (Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
          Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
          Kappa22*(Kappa20^2 + Kappa21^2 + Kappa22^2))*mDxbar222)) + 
    15*(Kappa00*(Kappa00*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
            Kappa02*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + 
            Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
          Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
            Kappa02*mDxbar222)) + Kappa10*
         (Kappa10*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
            Kappa02*mDxbar220) + Kappa11*(Kappa00*mDxbar201 + 
            Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
          Kappa12*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
            Kappa02*mDxbar222)) + Kappa20*
         (Kappa20*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
            Kappa02*mDxbar220) + Kappa21*(Kappa00*mDxbar201 + 
            Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
          Kappa22*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
            Kappa02*mDxbar222))) + 
      Kappa01*(Kappa01*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
            Kappa02*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + 
            Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
          Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
            Kappa02*mDxbar222)) + Kappa11*
         (Kappa10*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
            Kappa02*mDxbar220) + Kappa11*(Kappa00*mDxbar201 + 
            Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
          Kappa12*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
            Kappa02*mDxbar222)) + Kappa21*
         (Kappa20*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
            Kappa02*mDxbar220) + Kappa21*(Kappa00*mDxbar201 + 
            Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
          Kappa22*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
            Kappa02*mDxbar222))) + 
      Kappa02*(Kappa02*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
            Kappa02*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + 
            Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
          Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
            Kappa02*mDxbar222)) + Kappa12*
         (Kappa10*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
            Kappa02*mDxbar220) + Kappa11*(Kappa00*mDxbar201 + 
            Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
          Kappa12*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
            Kappa02*mDxbar222)) + Kappa22*
         (Kappa20*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
            Kappa02*mDxbar220) + Kappa21*(Kappa00*mDxbar201 + 
            Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
          Kappa22*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
            Kappa02*mDxbar222))) + 
      Kappa10*(Kappa00*(Kappa00*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
            Kappa12*mDxbar220) + Kappa01*(Kappa10*mDxbar201 + 
            Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
          Kappa02*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
            Kappa12*mDxbar222)) + Kappa10*
         (Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
            Kappa12*mDxbar220) + Kappa11*(Kappa10*mDxbar201 + 
            Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
          Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
            Kappa12*mDxbar222)) + Kappa20*
         (Kappa20*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
            Kappa12*mDxbar220) + Kappa21*(Kappa10*mDxbar201 + 
            Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
          Kappa22*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
            Kappa12*mDxbar222))) + 
      Kappa11*(Kappa01*(Kappa00*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
            Kappa12*mDxbar220) + Kappa01*(Kappa10*mDxbar201 + 
            Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
          Kappa02*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
            Kappa12*mDxbar222)) + Kappa11*
         (Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
            Kappa12*mDxbar220) + Kappa11*(Kappa10*mDxbar201 + 
            Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
          Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
            Kappa12*mDxbar222)) + Kappa21*
         (Kappa20*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
            Kappa12*mDxbar220) + Kappa21*(Kappa10*mDxbar201 + 
            Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
          Kappa22*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
            Kappa12*mDxbar222))) + 
      Kappa12*(Kappa02*(Kappa00*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
            Kappa12*mDxbar220) + Kappa01*(Kappa10*mDxbar201 + 
            Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
          Kappa02*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
            Kappa12*mDxbar222)) + Kappa12*
         (Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
            Kappa12*mDxbar220) + Kappa11*(Kappa10*mDxbar201 + 
            Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
          Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
            Kappa12*mDxbar222)) + Kappa22*
         (Kappa20*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
            Kappa12*mDxbar220) + Kappa21*(Kappa10*mDxbar201 + 
            Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
          Kappa22*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
            Kappa12*mDxbar222))) + 
      Kappa20*(Kappa00*(Kappa00*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
            Kappa22*mDxbar220) + Kappa01*(Kappa20*mDxbar201 + 
            Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
          Kappa02*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
            Kappa22*mDxbar222)) + Kappa10*
         (Kappa10*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
            Kappa22*mDxbar220) + Kappa11*(Kappa20*mDxbar201 + 
            Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
          Kappa12*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
            Kappa22*mDxbar222)) + Kappa20*
         (Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
            Kappa22*mDxbar220) + Kappa21*(Kappa20*mDxbar201 + 
            Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
          Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
            Kappa22*mDxbar222))) + 
      Kappa21*(Kappa01*(Kappa00*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
            Kappa22*mDxbar220) + Kappa01*(Kappa20*mDxbar201 + 
            Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
          Kappa02*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
            Kappa22*mDxbar222)) + Kappa11*
         (Kappa10*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
            Kappa22*mDxbar220) + Kappa11*(Kappa20*mDxbar201 + 
            Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
          Kappa12*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
            Kappa22*mDxbar222)) + Kappa21*
         (Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
            Kappa22*mDxbar220) + Kappa21*(Kappa20*mDxbar201 + 
            Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
          Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
            Kappa22*mDxbar222))) + 
      Kappa22*(Kappa02*(Kappa00*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
            Kappa22*mDxbar220) + Kappa01*(Kappa20*mDxbar201 + 
            Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
          Kappa02*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
            Kappa22*mDxbar222)) + Kappa12*
         (Kappa10*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
            Kappa22*mDxbar220) + Kappa11*(Kappa20*mDxbar201 + 
            Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
          Kappa12*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
            Kappa22*mDxbar222)) + Kappa22*
         (Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
            Kappa22*mDxbar220) + Kappa21*(Kappa20*mDxbar201 + 
            Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
          Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
            Kappa22*mDxbar222)))) - 
    3*g1^2*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
      Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + 
      Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
      Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211)) - 
    15*g2^2*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
      Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + 
      Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
      Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211)) + 
    20*(Lambda1200*(Lambda1200*(Lambda1200*(Lambda1200*mH1I200 + 
            Lambda1201*mH1I201) + Lambda1210*(Lambda1210*mH1I200 + 
            Lambda1211*mH1I201)) + Lambda1201*
         (Lambda1201*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
          Lambda1211*(Lambda1210*mH1I200 + Lambda1211*mH1I201))) + 
      Lambda1210*(Lambda1210*(Lambda1200*(Lambda1200*mH1I200 + 
            Lambda1201*mH1I201) + Lambda1210*(Lambda1210*mH1I200 + 
            Lambda1211*mH1I201)) + Lambda1211*
         (Lambda1201*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
          Lambda1211*(Lambda1210*mH1I200 + Lambda1211*mH1I201))) + 
      Lambda1201*(Lambda1200*(Lambda1200*(Lambda1200*mH1I210 + 
            Lambda1201*mH1I211) + Lambda1210*(Lambda1210*mH1I210 + 
            Lambda1211*mH1I211)) + Lambda1201*
         (Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
          Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211))) + 
      Lambda1211*(Lambda1210*(Lambda1200*(Lambda1200*mH1I210 + 
            Lambda1201*mH1I211) + Lambda1210*(Lambda1210*mH1I210 + 
            Lambda1211*mH1I211)) + Lambda1211*
         (Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
          Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211)))) - 
    3*g1^2*((Lambda1200^2 + Lambda1201^2)*mH2I200 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + 
      (Lambda1210^2 + Lambda1211^2)*mH2I211) - 
    15*g2^2*((Lambda1200^2 + Lambda1201^2)*mH2I200 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + 
      (Lambda1210^2 + Lambda1211^2)*mH2I211) + 
    10*((Lambda1200*(Lambda1200*(Lambda1200^2 + Lambda1201^2) + 
          Lambda1210*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)) + 
        Lambda1201*(Lambda1201*(Lambda1200^2 + Lambda1201^2) + 
          Lambda1211*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)))*
       mH2I200 + (Lambda1200*(Lambda1200*(Lambda1200*Lambda1210 + 
            Lambda1201*Lambda1211) + Lambda1210*(Lambda1210^2 + 
            Lambda1211^2)) + Lambda1201*(Lambda1201*(Lambda1200*Lambda1210 + 
            Lambda1201*Lambda1211) + Lambda1211*(Lambda1210^2 + 
            Lambda1211^2)))*mH2I201 + 
      (Lambda1210*(Lambda1200*(Lambda1200^2 + Lambda1201^2) + 
          Lambda1210*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)) + 
        Lambda1211*(Lambda1201*(Lambda1200^2 + Lambda1201^2) + 
          Lambda1211*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)))*
       mH2I210 + (Lambda1210*(Lambda1200*(Lambda1200*Lambda1210 + 
            Lambda1201*Lambda1211) + Lambda1210*(Lambda1210^2 + 
            Lambda1211^2)) + Lambda1211*(Lambda1201*(Lambda1200*Lambda1210 + 
            Lambda1201*Lambda1211) + Lambda1211*(Lambda1210^2 + 
            Lambda1211^2)))*mH2I211) + 
    10*(Lambda1200*(Lambda1200*((Lambda1200^2 + Lambda1201^2)*mH2I200 + 
          (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210) + 
        Lambda1210*((Lambda1200^2 + Lambda1201^2)*mH2I201 + 
          (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I211)) + 
      Lambda1201*(Lambda1201*((Lambda1200^2 + Lambda1201^2)*mH2I200 + 
          (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210) + 
        Lambda1211*((Lambda1200^2 + Lambda1201^2)*mH2I201 + 
          (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I211)) + 
      Lambda1210*(Lambda1200*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*
           mH2I200 + (Lambda1210^2 + Lambda1211^2)*mH2I210) + 
        Lambda1210*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
          (Lambda1210^2 + Lambda1211^2)*mH2I211)) + 
      Lambda1211*(Lambda1201*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*
           mH2I200 + (Lambda1210^2 + Lambda1211^2)*mH2I210) + 
        Lambda1211*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
          (Lambda1210^2 + Lambda1211^2)*mH2I211))) - 
    2*g1^2*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
      Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*ms2 - 
    40*g3^2*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
      Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*ms2 + 
    30*(Kappa00*(Kappa00*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
        Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
        Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
      Kappa01*(Kappa01*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
        Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
        Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
      Kappa02*(Kappa02*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
        Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
        Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
      Kappa10*(Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
          Kappa02*Kappa12) + Kappa10*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
        Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
      Kappa11*(Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
          Kappa02*Kappa12) + Kappa11*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
        Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
      Kappa12*(Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + 
          Kappa02*Kappa12) + Kappa12*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
        Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
      Kappa20*(Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
          Kappa02*Kappa22) + Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
          Kappa12*Kappa22) + Kappa20*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
      Kappa21*(Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
          Kappa02*Kappa22) + Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
          Kappa12*Kappa22) + Kappa21*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
      Kappa22*(Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 + 
          Kappa02*Kappa22) + Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + 
          Kappa12*Kappa22) + Kappa22*(Kappa20^2 + Kappa21^2 + Kappa22^2)))*
     ms2 - 3*g1^2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*
     ms2 - 15*g2^2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
      Lambda1211^2)*ms2 + 
    20*(Lambda1200*(Lambda1200*(Lambda1200^2 + Lambda1201^2) + 
        Lambda1210*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)) + 
      Lambda1201*(Lambda1201*(Lambda1200^2 + Lambda1201^2) + 
        Lambda1211*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)) + 
      Lambda1210*(Lambda1200*(Lambda1200*Lambda1210 + 
          Lambda1201*Lambda1211) + Lambda1210*(Lambda1210^2 + 
          Lambda1211^2)) + Lambda1211*(Lambda1201*(Lambda1200*Lambda1210 + 
          Lambda1201*Lambda1211) + Lambda1211*(Lambda1210^2 + Lambda1211^2)))*
     ms2 + 20*Lambdax^4*(mHd2 + mHu2 + ms2) - 
    15*gN^2*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx202 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
      (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx212 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 + 
      (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222)*QDxbarp^2 - 
    15*gN^2*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
        Kappa02*mDxbar220) + Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
        Kappa12*mDxbar220) + Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
        Kappa22*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + 
        Kappa02*mDxbar221) + Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + 
        Kappa12*mDxbar221) + Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + 
        Kappa22*mDxbar221) + Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
        Kappa02*mDxbar222) + Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
        Kappa12*mDxbar222) + Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
        Kappa22*mDxbar222))*QDxbarp^2 - 
    15*gN^2*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
      Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*ms2*QDxbarp^2 - 
    15*gN^2*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx202 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
      (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx212 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 + 
      (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222)*QDxp^2 - 
    15*gN^2*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
        Kappa02*mDxbar220) + Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
        Kappa12*mDxbar220) + Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
        Kappa22*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + 
        Kappa02*mDxbar221) + Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + 
        Kappa12*mDxbar221) + Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + 
        Kappa22*mDxbar221) + Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
        Kappa02*mDxbar222) + Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
        Kappa12*mDxbar222) + Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
        Kappa22*mDxbar222))*QDxp^2 - 15*gN^2*(Kappa00^2 + Kappa01^2 + 
      Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + 
      Kappa22^2)*ms2*QDxp^2 - 
    10*gN^2*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
      Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + 
      Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
      Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211))*QH1p^2 - 
    10*gN^2*((Lambda1200^2 + Lambda1201^2)*mH2I200 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + 
      (Lambda1210^2 + Lambda1211^2)*mH2I211)*QH1p^2 - 
    10*gN^2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*ms2*
     QH1p^2 - 10*gN^2*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
      Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + 
      Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
      Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211))*QH2p^2 - 
    10*gN^2*((Lambda1200^2 + Lambda1201^2)*mH2I200 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + 
      (Lambda1210^2 + Lambda1211^2)*mH2I211)*QH2p^2 - 
    10*gN^2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*ms2*
     QH2p^2 + 15*gN^2*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx202 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
      (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx212 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 + 
      (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222)*QSp^2 + 
    15*gN^2*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
        Kappa02*mDxbar220) + Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
        Kappa12*mDxbar220) + Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
        Kappa22*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + 
        Kappa02*mDxbar221) + Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + 
        Kappa12*mDxbar221) + Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + 
        Kappa22*mDxbar221) + Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
        Kappa02*mDxbar222) + Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
        Kappa12*mDxbar222) + Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
        Kappa22*mDxbar222))*QSp^2 + 
    10*gN^2*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
      Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + 
      Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
      Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211))*QSp^2 + 
    10*gN^2*((Lambda1200^2 + Lambda1201^2)*mH2I200 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + 
      (Lambda1210^2 + Lambda1211^2)*mH2I211)*QSp^2 + 
    15*gN^2*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
      Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*ms2*QSp^2 + 
    10*gN^2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*ms2*
     QSp^2 - 10*gN^4*QSp^2*(3*(md200 + md211 + md222)*Qdp^2 + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp^2 + 
      3*(mDx200 + mDx211 + mDx222)*QDxp^2 + (me200 + me211 + me222)*Qep^2 + 
      2*(mH1I200 + mH1I211)*QH1p^2 + 2*mHd2*QH1p^2 + 
      2*(mH2I200 + mH2I211)*QH2p^2 + 2*mHu2*QH2p^2 + 2*mHpbar2*QHpbarp^2 + 
      2*mHp2*QHpp^2 + 2*(ml200 + ml211 + ml222)*QLp^2 + 
      6*(mq200 + mq211 + mq222)*QQp^2 + ms2*QSp^2 + (msI200 + msI211)*QSp^2 + 
      3*(mu200 + mu211 + mu222)*Qup^2) + 4*g1^2*MassB*
     (Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
      Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
      Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
    80*g3^2*MassG*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
      Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
      Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
    15*gN^2*MassBp*QDxbarp^2*(Kappa00*TKappa00 + Kappa01*TKappa01 + 
      Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + 
      Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + 
      Kappa22*TKappa22) + 15*gN^2*MassBp*QDxp^2*(Kappa00*TKappa00 + 
      Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + 
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + 
      Kappa21*TKappa21 + Kappa22*TKappa22) - 15*gN^2*MassBp*QSp^2*
     (Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
      Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
      Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) - 
    2*g1^2*(TKappa00^2 + TKappa01^2 + TKappa02^2 + TKappa10^2 + TKappa11^2 + 
      TKappa12^2 + TKappa20^2 + TKappa21^2 + TKappa22^2) - 
    40*g3^2*(TKappa00^2 + TKappa01^2 + TKappa02^2 + TKappa10^2 + TKappa11^2 + 
      TKappa12^2 + TKappa20^2 + TKappa21^2 + TKappa22^2) - 
    15*gN^2*QDxbarp^2*(TKappa00^2 + TKappa01^2 + TKappa02^2 + TKappa10^2 + 
      TKappa11^2 + TKappa12^2 + TKappa20^2 + TKappa21^2 + TKappa22^2) - 
    15*gN^2*QDxp^2*(TKappa00^2 + TKappa01^2 + TKappa02^2 + TKappa10^2 + 
      TKappa11^2 + TKappa12^2 + TKappa20^2 + TKappa21^2 + TKappa22^2) + 
    15*gN^2*QSp^2*(TKappa00^2 + TKappa01^2 + TKappa02^2 + TKappa10^2 + 
      TKappa11^2 + TKappa12^2 + TKappa20^2 + TKappa21^2 + TKappa22^2) + 
    30*(TKappa00*((Kappa00^2 + Kappa01^2 + Kappa02^2)*TKappa00 + 
        (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa10 + 
        (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa20) + 
      TKappa10*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
         TKappa00 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*TKappa10 + 
        (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa20) + 
      TKappa20*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
         TKappa00 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
         TKappa10 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa20) + 
      TKappa01*((Kappa00^2 + Kappa01^2 + Kappa02^2)*TKappa01 + 
        (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa11 + 
        (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa21) + 
      TKappa11*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
         TKappa01 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*TKappa11 + 
        (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa21) + 
      TKappa21*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
         TKappa01 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
         TKappa11 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa21) + 
      TKappa02*((Kappa00^2 + Kappa01^2 + Kappa02^2)*TKappa02 + 
        (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa12 + 
        (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa22) + 
      TKappa12*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
         TKappa02 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*TKappa12 + 
        (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa22) + 
      TKappa22*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
         TKappa02 + (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
         TKappa12 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa22)) + 
    30*(Kappa00*(TKappa00*(Kappa00*TKappa00 + Kappa01*TKappa01 + 
          Kappa02*TKappa02) + TKappa10*(Kappa00*TKappa10 + Kappa01*TKappa11 + 
          Kappa02*TKappa12) + TKappa20*(Kappa00*TKappa20 + Kappa01*TKappa21 + 
          Kappa02*TKappa22)) + Kappa01*
       (TKappa01*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02) + 
        TKappa11*(Kappa00*TKappa10 + Kappa01*TKappa11 + Kappa02*TKappa12) + 
        TKappa21*(Kappa00*TKappa20 + Kappa01*TKappa21 + Kappa02*TKappa22)) + 
      Kappa02*(TKappa02*(Kappa00*TKappa00 + Kappa01*TKappa01 + 
          Kappa02*TKappa02) + TKappa12*(Kappa00*TKappa10 + Kappa01*TKappa11 + 
          Kappa02*TKappa12) + TKappa22*(Kappa00*TKappa20 + Kappa01*TKappa21 + 
          Kappa02*TKappa22)) + Kappa10*
       (TKappa00*(Kappa10*TKappa00 + Kappa11*TKappa01 + Kappa12*TKappa02) + 
        TKappa10*(Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12) + 
        TKappa20*(Kappa10*TKappa20 + Kappa11*TKappa21 + Kappa12*TKappa22)) + 
      Kappa11*(TKappa01*(Kappa10*TKappa00 + Kappa11*TKappa01 + 
          Kappa12*TKappa02) + TKappa11*(Kappa10*TKappa10 + Kappa11*TKappa11 + 
          Kappa12*TKappa12) + TKappa21*(Kappa10*TKappa20 + Kappa11*TKappa21 + 
          Kappa12*TKappa22)) + Kappa12*
       (TKappa02*(Kappa10*TKappa00 + Kappa11*TKappa01 + Kappa12*TKappa02) + 
        TKappa12*(Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12) + 
        TKappa22*(Kappa10*TKappa20 + Kappa11*TKappa21 + Kappa12*TKappa22)) + 
      Kappa20*(TKappa00*(Kappa20*TKappa00 + Kappa21*TKappa01 + 
          Kappa22*TKappa02) + TKappa10*(Kappa20*TKappa10 + Kappa21*TKappa11 + 
          Kappa22*TKappa12) + TKappa20*(Kappa20*TKappa20 + Kappa21*TKappa21 + 
          Kappa22*TKappa22)) + Kappa21*
       (TKappa01*(Kappa20*TKappa00 + Kappa21*TKappa01 + Kappa22*TKappa02) + 
        TKappa11*(Kappa20*TKappa10 + Kappa21*TKappa11 + Kappa22*TKappa12) + 
        TKappa21*(Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)) + 
      Kappa22*(TKappa02*(Kappa20*TKappa00 + Kappa21*TKappa01 + 
          Kappa22*TKappa02) + TKappa12*(Kappa20*TKappa10 + Kappa21*TKappa11 + 
          Kappa22*TKappa12) + TKappa22*(Kappa20*TKappa20 + Kappa21*TKappa21 + 
          Kappa22*TKappa22))) + 6*g1^2*MassB*(Lambda1200*TLambda1200 + 
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + 
      Lambda1211*TLambda1211) + 30*g2^2*MassWB*(Lambda1200*TLambda1200 + 
      Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + 
      Lambda1211*TLambda1211) + 10*gN^2*MassBp*QH1p^2*
     (Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
    10*gN^2*MassBp*QH2p^2*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211) - 
    10*gN^2*MassBp*QSp^2*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211) - 
    3*g1^2*(TLambda1200^2 + TLambda1201^2 + TLambda1210^2 + TLambda1211^2) - 
    15*g2^2*(TLambda1200^2 + TLambda1201^2 + TLambda1210^2 + TLambda1211^2) - 
    10*gN^2*QH1p^2*(TLambda1200^2 + TLambda1201^2 + TLambda1210^2 + 
      TLambda1211^2) - 10*gN^2*QH2p^2*(TLambda1200^2 + TLambda1201^2 + 
      TLambda1210^2 + TLambda1211^2) + 10*gN^2*QSp^2*
     (TLambda1200^2 + TLambda1201^2 + TLambda1210^2 + TLambda1211^2) + 
    20*(TLambda1200*((Lambda1200^2 + Lambda1201^2)*TLambda1200 + 
        (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1210) + 
      TLambda1210*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*
         TLambda1200 + (Lambda1210^2 + Lambda1211^2)*TLambda1210) + 
      TLambda1201*((Lambda1200^2 + Lambda1201^2)*TLambda1201 + 
        (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1211) + 
      TLambda1211*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*
         TLambda1201 + (Lambda1210^2 + Lambda1211^2)*TLambda1211)) + 
    20*(Lambda1200*(TLambda1200*(Lambda1200*TLambda1200 + 
          Lambda1201*TLambda1201) + TLambda1210*(Lambda1200*TLambda1210 + 
          Lambda1201*TLambda1211)) + Lambda1201*
       (TLambda1201*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201) + 
        TLambda1211*(Lambda1200*TLambda1210 + Lambda1201*TLambda1211)) + 
      Lambda1210*(TLambda1200*(Lambda1210*TLambda1200 + 
          Lambda1211*TLambda1201) + TLambda1210*(Lambda1210*TLambda1210 + 
          Lambda1211*TLambda1211)) + Lambda1211*
       (TLambda1201*(Lambda1210*TLambda1200 + Lambda1211*TLambda1201) + 
        TLambda1211*(Lambda1210*TLambda1210 + Lambda1211*TLambda1211))) - 
    5*gN^2*MassBp*(4*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
        Lambda1211^2)*MassBp*QH1p^2 + 4*(Lambda1200^2 + Lambda1201^2 + 
        Lambda1210^2 + Lambda1211^2)*MassBp*QH2p^2 - 
      4*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*MassBp*
       QSp^2 + 54*gN^2*MassBp*Qdp^2*QSp^2 + 54*gN^2*MassBp*QDxbarp^2*QSp^2 + 
      54*gN^2*MassBp*QDxp^2*QSp^2 + 18*gN^2*MassBp*Qep^2*QSp^2 + 
      36*gN^2*MassBp*QH1p^2*QSp^2 + 36*gN^2*MassBp*QH2p^2*QSp^2 + 
      12*gN^2*MassBp*QHpbarp^2*QSp^2 + 12*gN^2*MassBp*QHpp^2*QSp^2 + 
      36*gN^2*MassBp*QLp^2*QSp^2 + 108*gN^2*MassBp*QQp^2*QSp^2 + 
      30*gN^2*MassBp*QSp^4 + 6*(Kappa00^2 + Kappa01^2 + Kappa02^2 + 
        Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + 
        Kappa22^2)*MassBp*(QDxbarp^2 + QDxp^2 - QSp^2) + 
      54*gN^2*MassBp*QSp^2*Qup^2 - 3*QDxbarp^2*(Kappa00*TKappa00 + 
        Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + 
        Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + 
        Kappa21*TKappa21 + Kappa22*TKappa22) - 
      3*QDxp^2*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      3*QSp^2*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) - 
      2*QH1p^2*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) - 
      2*QH2p^2*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      2*QSp^2*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      2*Lambdax*(QH1p^2 + QH2p^2 - QSp^2)*(2*Lambdax*MassBp - TLambdax)) + 
    3*g1^2*Lambdax*MassB*TLambdax + 15*g2^2*Lambdax*MassWB*TLambdax + 
    10*gN^2*Lambdax*MassBp*QH1p^2*TLambdax + 10*gN^2*Lambdax*MassBp*QH2p^2*
     TLambdax - 10*gN^2*Lambdax*MassBp*QSp^2*TLambdax - 3*g1^2*TLambdax^2 - 
    15*g2^2*TLambdax^2 - 10*gN^2*QH1p^2*TLambdax^2 - 
    10*gN^2*QH2p^2*TLambdax^2 + 10*gN^2*QSp^2*TLambdax^2 + 
    15*Lambdax*TLambdax*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
      TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
    15*TLambdax^2*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
      Yd20^2 + Yd21^2 + Yd22^2) + 5*Lambdax*TLambdax*
     (TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
      TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 
    5*TLambdax^2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + 
      Ye20^2 + Ye21^2 + Ye22^2) + 15*Lambdax*TLambdax*
     (TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + TYu11*Yu11 + 
      TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 
    15*TLambdax^2*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + 
      Yu20^2 + Yu21^2 + Yu22^2) + Lambdax*(-3*g1^2*Lambdax*mHd2 - 
      15*g2^2*Lambdax*mHd2 - 3*g1^2*Lambdax*mHu2 - 15*g2^2*Lambdax*mHu2 - 
      3*g1^2*Lambdax*ms2 - 15*g2^2*Lambdax*ms2 - 10*gN^2*Lambdax*mHd2*
       QH1p^2 - 10*gN^2*Lambdax*mHu2*QH1p^2 - 10*gN^2*Lambdax*ms2*QH1p^2 - 
      10*gN^2*Lambdax*mHd2*QH2p^2 - 10*gN^2*Lambdax*mHu2*QH2p^2 - 
      10*gN^2*Lambdax*ms2*QH2p^2 + 10*gN^2*Lambdax*mHd2*QSp^2 + 
      10*gN^2*Lambdax*mHu2*QSp^2 + 10*gN^2*Lambdax*ms2*QSp^2 + 
      40*Lambdax*TLambdax^2 + 3*g1^2*MassB*(-2*Lambdax*MassB + TLambdax) + 
      15*g2^2*MassWB*(-2*Lambdax*MassWB + TLambdax) + 
      15*Lambdax*(TYd00^2 + TYd01^2 + TYd02^2 + TYd10^2 + TYd11^2 + TYd12^2 + 
        TYd20^2 + TYd21^2 + TYd22^2) + 5*Lambdax*(TYe00^2 + TYe01^2 + 
        TYe02^2 + TYe10^2 + TYe11^2 + TYe12^2 + TYe20^2 + TYe21^2 + 
        TYe22^2) + 15*Lambdax*(TYu00^2 + TYu01^2 + TYu02^2 + TYu10^2 + 
        TYu11^2 + TYu12^2 + TYu20^2 + TYu21^2 + TYu22^2) + 
      15*TLambdax*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
        TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      30*Lambdax*mHd2*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
        Yd20^2 + Yd21^2 + Yd22^2) + 15*Lambdax*mHu2*(Yd00^2 + Yd01^2 + 
        Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
      15*Lambdax*ms2*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
        Yd20^2 + Yd21^2 + Yd22^2) + 15*Lambdax*
       (Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
        Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
        Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
        Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
        Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
        Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
        Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22) + 
        Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22) + 
        Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
      15*Lambdax*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + 
        Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + 
        Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) + 
      5*TLambdax*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
        TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 
      10*Lambdax*mHd2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + 
        Ye20^2 + Ye21^2 + Ye22^2) + 5*Lambdax*mHu2*(Ye00^2 + Ye01^2 + 
        Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
      5*Lambdax*ms2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + 
        Ye20^2 + Ye21^2 + Ye22^2) + 5*Lambdax*
       (Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
        Ye10*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
        Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
        Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
        Ye11*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
        Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
        Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22) + 
        Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22) + 
        Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
      5*Lambdax*(Ye00*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
        Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
        Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
        Ye10*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
        Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
        Ye12*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
        Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22) + 
        Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + 
        Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)) + 
      15*TLambdax*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + 
        TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 
      15*Lambdax*mHd2*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + 
        Yu20^2 + Yu21^2 + Yu22^2) + 30*Lambdax*mHu2*(Yu00^2 + Yu01^2 + 
        Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + Yu20^2 + Yu21^2 + Yu22^2) + 
      15*Lambdax*ms2*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + 
        Yu20^2 + Yu21^2 + Yu22^2) + 15*Lambdax*
       (Yu00*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
        Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
        Yu02*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
        Yu10*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
        Yu11*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
        Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
        Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22) + 
        Yu21*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22) + 
        Yu22*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) + 
      15*Lambdax*(Yu00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
        Yu10*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
        Yu20*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
        Yu01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
        Yu11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
        Yu21*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
        Yu02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22) + 
        Yu12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22) + 
        Yu22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22))) - 
    gN^2*QSp*(2*(md200 + md211 + md222)*Qdp*(g1^2 + 20*g3^2 + 
        15*gN^2*Qdp^2) + 2*g1^2*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      40*g3^2*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp - 
      15*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar201 + 
          Kappa02*mDxbar202) + Kappa10*(Kappa10*mDxbar200 + 
          Kappa11*mDxbar201 + Kappa12*mDxbar202) + 
        Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar201 + Kappa22*mDxbar202) + 
        Kappa01*(Kappa00*mDxbar210 + Kappa01*mDxbar211 + Kappa02*mDxbar212) + 
        Kappa11*(Kappa10*mDxbar210 + Kappa11*mDxbar211 + Kappa12*mDxbar212) + 
        Kappa21*(Kappa20*mDxbar210 + Kappa21*mDxbar211 + Kappa22*mDxbar212) + 
        Kappa02*(Kappa00*mDxbar220 + Kappa01*mDxbar221 + Kappa02*mDxbar222) + 
        Kappa12*(Kappa10*mDxbar220 + Kappa11*mDxbar221 + Kappa12*mDxbar222) + 
        Kappa22*(Kappa20*mDxbar220 + Kappa21*mDxbar221 + Kappa22*mDxbar222))*
       QDxbarp + 30*gN^2*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp^3 + 
      2*g1^2*(mDx200 + mDx211 + mDx222)*QDxp + 
      40*g3^2*(mDx200 + mDx211 + mDx222)*QDxp - 
      15*(Kappa00*(Kappa00*mDx200 + Kappa10*mDx201 + Kappa20*mDx202) + 
        Kappa01*(Kappa01*mDx200 + Kappa11*mDx201 + Kappa21*mDx202) + 
        Kappa02*(Kappa02*mDx200 + Kappa12*mDx201 + Kappa22*mDx202) + 
        Kappa10*(Kappa00*mDx210 + Kappa10*mDx211 + Kappa20*mDx212) + 
        Kappa11*(Kappa01*mDx210 + Kappa11*mDx211 + Kappa21*mDx212) + 
        Kappa12*(Kappa02*mDx210 + Kappa12*mDx211 + Kappa22*mDx212) + 
        Kappa20*(Kappa00*mDx220 + Kappa10*mDx221 + Kappa20*mDx222) + 
        Kappa21*(Kappa01*mDx220 + Kappa11*mDx221 + Kappa21*mDx222) + 
        Kappa22*(Kappa02*mDx220 + Kappa12*mDx221 + Kappa22*mDx222))*QDxp + 
      30*gN^2*(mDx200 + mDx211 + mDx222)*QDxp^3 + 
      6*g1^2*(me200 + me211 + me222)*Qep + 10*gN^2*(me200 + me211 + me222)*
       Qep^3 + 3*g1^2*(mH1I200 + mH1I211)*QH1p + 15*g2^2*(mH1I200 + mH1I211)*
       QH1p - 10*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I210) + 
        Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I210) + 
        Lambda1201*(Lambda1200*mH1I201 + Lambda1201*mH1I211) + 
        Lambda1211*(Lambda1210*mH1I201 + Lambda1211*mH1I211))*QH1p + 
      3*g1^2*mHd2*QH1p + 15*g2^2*mHd2*QH1p + 20*gN^2*(mH1I200 + mH1I211)*
       QH1p^3 + 20*gN^2*mHd2*QH1p^3 + 3*g1^2*(mH2I200 + mH2I211)*QH2p + 
      15*g2^2*(mH2I200 + mH2I211)*QH2p - 
      10*(Lambda1200*(Lambda1200*mH2I200 + Lambda1210*mH2I201) + 
        Lambda1201*(Lambda1201*mH2I200 + Lambda1211*mH2I201) + 
        Lambda1210*(Lambda1200*mH2I210 + Lambda1210*mH2I211) + 
        Lambda1211*(Lambda1201*mH2I210 + Lambda1211*mH2I211))*QH2p + 
      3*g1^2*mHu2*QH2p + 15*g2^2*mHu2*QH2p + 20*gN^2*(mH2I200 + mH2I211)*
       QH2p^3 + 20*gN^2*mHu2*QH2p^3 + 3*g1^2*mHpbar2*QHpbarp + 
      15*g2^2*mHpbar2*QHpbarp + 20*gN^2*mHpbar2*QHpbarp^3 + 
      3*g1^2*mHp2*QHpp + 15*g2^2*mHp2*QHpp + 20*gN^2*mHp2*QHpp^3 + 
      3*g1^2*(ml200 + ml211 + ml222)*QLp + 15*g2^2*(ml200 + ml211 + ml222)*
       QLp + 20*gN^2*(ml200 + ml211 + ml222)*QLp^3 + 
      g1^2*(mq200 + mq211 + mq222)*QQp + 45*g2^2*(mq200 + mq211 + mq222)*
       QQp + 80*g3^2*(mq200 + mq211 + mq222)*QQp + 
      60*gN^2*(mq200 + mq211 + mq222)*QQp^3 - 
      15*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
        Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*ms2*QSp - 
      10*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*ms2*
       QSp + 10*gN^2*ms2*QSp^3 + 10*gN^2*(msI200 + msI211)*QSp^3 - 
      10*Lambdax^2*(mHd2*QH1p + mHu2*QH2p + ms2*QSp) + 
      8*g1^2*(mu200 + mu211 + mu222)*Qup + 40*g3^2*(mu200 + mu211 + mu222)*
       Qup + 30*gN^2*(mu200 + mu211 + mu222)*Qup^3 - 
      30*mHd2*QH1p*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
        Yd20^2 + Yd21^2 + Yd22^2) - 
      30*QQp*(Yd00*(mq200*Yd00 + mq210*Yd01 + mq220*Yd02) + 
        Yd01*(mq201*Yd00 + mq211*Yd01 + mq221*Yd02) + 
        Yd02*(mq202*Yd00 + mq212*Yd01 + mq222*Yd02) + 
        Yd10*(mq200*Yd10 + mq210*Yd11 + mq220*Yd12) + 
        Yd11*(mq201*Yd10 + mq211*Yd11 + mq221*Yd12) + 
        Yd12*(mq202*Yd10 + mq212*Yd11 + mq222*Yd12) + 
        Yd20*(mq200*Yd20 + mq210*Yd21 + mq220*Yd22) + 
        Yd21*(mq201*Yd20 + mq211*Yd21 + mq221*Yd22) + 
        Yd22*(mq202*Yd20 + mq212*Yd21 + mq222*Yd22)) - 
      30*Qdp*(md200*(Yd00^2 + Yd01^2 + Yd02^2) + 
        md201*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
        md210*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
        md211*(Yd10^2 + Yd11^2 + Yd12^2) + md202*(Yd00*Yd20 + Yd01*Yd21 + 
          Yd02*Yd22) + md220*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
        md212*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
        md221*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
        md222*(Yd20^2 + Yd21^2 + Yd22^2)) - 10*mHd2*QH1p*
       (Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
        Ye21^2 + Ye22^2) - 10*QLp*(Ye00*(ml200*Ye00 + ml210*Ye01 + 
          ml220*Ye02) + Ye01*(ml201*Ye00 + ml211*Ye01 + ml221*Ye02) + 
        Ye02*(ml202*Ye00 + ml212*Ye01 + ml222*Ye02) + 
        Ye10*(ml200*Ye10 + ml210*Ye11 + ml220*Ye12) + 
        Ye11*(ml201*Ye10 + ml211*Ye11 + ml221*Ye12) + 
        Ye12*(ml202*Ye10 + ml212*Ye11 + ml222*Ye12) + 
        Ye20*(ml200*Ye20 + ml210*Ye21 + ml220*Ye22) + 
        Ye21*(ml201*Ye20 + ml211*Ye21 + ml221*Ye22) + 
        Ye22*(ml202*Ye20 + ml212*Ye21 + ml222*Ye22)) - 
      10*Qep*(me200*(Ye00^2 + Ye01^2 + Ye02^2) + 
        me201*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
        me210*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
        me211*(Ye10^2 + Ye11^2 + Ye12^2) + me202*(Ye00*Ye20 + Ye01*Ye21 + 
          Ye02*Ye22) + me220*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
        me212*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
        me221*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
        me222*(Ye20^2 + Ye21^2 + Ye22^2)) - 30*mHu2*QH2p*
       (Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + Yu20^2 + 
        Yu21^2 + Yu22^2) - 30*QQp*(Yu00*(mq200*Yu00 + mq210*Yu01 + 
          mq220*Yu02) + Yu01*(mq201*Yu00 + mq211*Yu01 + mq221*Yu02) + 
        Yu02*(mq202*Yu00 + mq212*Yu01 + mq222*Yu02) + 
        Yu10*(mq200*Yu10 + mq210*Yu11 + mq220*Yu12) + 
        Yu11*(mq201*Yu10 + mq211*Yu11 + mq221*Yu12) + 
        Yu12*(mq202*Yu10 + mq212*Yu11 + mq222*Yu12) + 
        Yu20*(mq200*Yu20 + mq210*Yu21 + mq220*Yu22) + 
        Yu21*(mq201*Yu20 + mq211*Yu21 + mq221*Yu22) + 
        Yu22*(mq202*Yu20 + mq212*Yu21 + mq222*Yu22)) - 
      30*Qup*(mu200*(Yu00^2 + Yu01^2 + Yu02^2) + 
        mu201*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + 
        mu210*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + 
        mu211*(Yu10^2 + Yu11^2 + Yu12^2) + mu202*(Yu00*Yu20 + Yu01*Yu21 + 
          Yu02*Yu22) + mu220*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 
        mu212*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 
        mu221*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 
        mu222*(Yu20^2 + Yu21^2 + Yu22^2)))))/5, 
 (6*(2*Kappa00*mDx200 + Kappa10*mDx201 + Kappa20*mDx202 + Kappa10*mDx210 + 
      Kappa20*mDx220) + 6*(2*Kappa00*mDxbar200 + Kappa01*mDxbar201 + 
      Kappa02*mDxbar202 + Kappa01*mDxbar210 + Kappa02*mDxbar220) + 
    12*Kappa00*ms2)*(2*(Kappa00*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
      Kappa10*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
      Kappa20*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
    Kappa00*((-4*g1^2)/15 - (16*g3^2)/3 + 3*(Kappa00^2 + Kappa01^2 + 
        Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + 
        Kappa21^2 + Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + 
        Lambda1210^2 + Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QDxbarp^2 - 
      2*gN^2*QDxp^2 - 2*gN^2*QSp^2)) + 
  (6*(2*Kappa01*mDx200 + Kappa11*mDx201 + Kappa21*mDx202 + Kappa11*mDx210 + 
      Kappa21*mDx220) + 6*(Kappa00*mDxbar201 + Kappa00*mDxbar210 + 
      2*Kappa01*mDxbar211 + Kappa02*mDxbar212 + Kappa02*mDxbar221) + 
    12*Kappa01*ms2)*(2*(Kappa01*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
      Kappa11*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
      Kappa21*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
    Kappa01*((-4*g1^2)/15 - (16*g3^2)/3 + 3*(Kappa00^2 + Kappa01^2 + 
        Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + 
        Kappa21^2 + Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + 
        Lambda1210^2 + Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QDxbarp^2 - 
      2*gN^2*QDxp^2 - 2*gN^2*QSp^2)) + 
  (6*(2*Kappa02*mDx200 + Kappa12*mDx201 + Kappa22*mDx202 + Kappa12*mDx210 + 
      Kappa22*mDx220) + 6*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
      Kappa00*mDxbar220 + Kappa01*mDxbar221 + 2*Kappa02*mDxbar222) + 
    12*Kappa02*ms2)*(2*(Kappa02*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 
      Kappa12*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
      Kappa22*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)) + 
    Kappa02*((-4*g1^2)/15 - (16*g3^2)/3 + 3*(Kappa00^2 + Kappa01^2 + 
        Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + 
        Kappa21^2 + Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + 
        Lambda1210^2 + Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QDxbarp^2 - 
      2*gN^2*QDxp^2 - 2*gN^2*QSp^2)) + 
  (6*(Kappa00*mDx201 + Kappa00*mDx210 + 2*Kappa10*mDx211 + Kappa20*mDx212 + 
      Kappa20*mDx221) + 6*(2*Kappa10*mDxbar200 + Kappa11*mDxbar201 + 
      Kappa12*mDxbar202 + Kappa11*mDxbar210 + Kappa12*mDxbar220) + 
    12*Kappa10*ms2)*
   (2*(Kappa00*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
      Kappa10*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
      Kappa20*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
    Kappa10*((-4*g1^2)/15 - (16*g3^2)/3 + 3*(Kappa00^2 + Kappa01^2 + 
        Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + 
        Kappa21^2 + Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + 
        Lambda1210^2 + Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QDxbarp^2 - 
      2*gN^2*QDxp^2 - 2*gN^2*QSp^2)) + 
  (6*(Kappa01*mDx201 + Kappa01*mDx210 + 2*Kappa11*mDx211 + Kappa21*mDx212 + 
      Kappa21*mDx221) + 6*(Kappa10*mDxbar201 + Kappa10*mDxbar210 + 
      2*Kappa11*mDxbar211 + Kappa12*mDxbar212 + Kappa12*mDxbar221) + 
    12*Kappa11*ms2)*
   (2*(Kappa01*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
      Kappa11*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
      Kappa21*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
    Kappa11*((-4*g1^2)/15 - (16*g3^2)/3 + 3*(Kappa00^2 + Kappa01^2 + 
        Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + 
        Kappa21^2 + Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + 
        Lambda1210^2 + Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QDxbarp^2 - 
      2*gN^2*QDxp^2 - 2*gN^2*QSp^2)) + 
  (6*(Kappa02*mDx201 + Kappa02*mDx210 + 2*Kappa12*mDx211 + Kappa22*mDx212 + 
      Kappa22*mDx221) + 6*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
      Kappa10*mDxbar220 + Kappa11*mDxbar221 + 2*Kappa12*mDxbar222) + 
    12*Kappa12*ms2)*
   (2*(Kappa02*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12) + 
      Kappa12*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 
      Kappa22*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)) + 
    Kappa12*((-4*g1^2)/15 - (16*g3^2)/3 + 3*(Kappa00^2 + Kappa01^2 + 
        Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + 
        Kappa21^2 + Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + 
        Lambda1210^2 + Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QDxbarp^2 - 
      2*gN^2*QDxp^2 - 2*gN^2*QSp^2)) + 
  (6*(Kappa00*mDx202 + Kappa10*mDx212 + Kappa00*mDx220 + Kappa10*mDx221 + 
      2*Kappa20*mDx222) + 6*(2*Kappa20*mDxbar200 + Kappa21*mDxbar201 + 
      Kappa22*mDxbar202 + Kappa21*mDxbar210 + Kappa22*mDxbar220) + 
    12*Kappa20*ms2)*
   (2*(Kappa00*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
      Kappa10*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
      Kappa20*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
    Kappa20*((-4*g1^2)/15 - (16*g3^2)/3 + 3*(Kappa00^2 + Kappa01^2 + 
        Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + 
        Kappa21^2 + Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + 
        Lambda1210^2 + Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QDxbarp^2 - 
      2*gN^2*QDxp^2 - 2*gN^2*QSp^2)) + 
  (6*(Kappa01*mDx202 + Kappa11*mDx212 + Kappa01*mDx220 + Kappa11*mDx221 + 
      2*Kappa21*mDx222) + 6*(Kappa20*mDxbar201 + Kappa20*mDxbar210 + 
      2*Kappa21*mDxbar211 + Kappa22*mDxbar212 + Kappa22*mDxbar221) + 
    12*Kappa21*ms2)*
   (2*(Kappa01*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
      Kappa11*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
      Kappa21*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
    Kappa21*((-4*g1^2)/15 - (16*g3^2)/3 + 3*(Kappa00^2 + Kappa01^2 + 
        Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + 
        Kappa21^2 + Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + 
        Lambda1210^2 + Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QDxbarp^2 - 
      2*gN^2*QDxp^2 - 2*gN^2*QSp^2)) + 
  (6*(Kappa02*mDx202 + Kappa12*mDx212 + Kappa02*mDx220 + Kappa12*mDx221 + 
      2*Kappa22*mDx222) + 6*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
      Kappa20*mDxbar220 + Kappa21*mDxbar221 + 2*Kappa22*mDxbar222) + 
    12*Kappa22*ms2)*
   (2*(Kappa02*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22) + 
      Kappa12*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22) + 
      Kappa22*(Kappa20^2 + Kappa21^2 + Kappa22^2)) + 
    Kappa22*((-4*g1^2)/15 - (16*g3^2)/3 + 3*(Kappa00^2 + Kappa01^2 + 
        Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + 
        Kappa21^2 + Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + 
        Lambda1210^2 + Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QDxbarp^2 - 
      2*gN^2*QDxp^2 - 2*gN^2*QSp^2)) + 
  (4*(2*Lambda1200*mH1I200 + Lambda1201*mH1I201 + Lambda1201*mH1I210) + 
    4*(2*Lambda1200*mH2I200 + Lambda1210*mH2I201 + Lambda1210*mH2I210) + 
    8*Lambda1200*ms2)*(2*(Lambda1200*(Lambda1200^2 + Lambda1201^2) + 
      Lambda1210*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)) + 
    Lambda1200*((-3*g1^2)/5 - 3*g2^2 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + 
        Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + 
        Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
        Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QH1p^2 - 2*gN^2*QH2p^2 - 
      2*gN^2*QSp^2)) + 
  (4*(Lambda1200*mH1I201 + Lambda1200*mH1I210 + 2*Lambda1201*mH1I211) + 
    4*(2*Lambda1201*mH2I200 + Lambda1211*mH2I201 + Lambda1211*mH2I210) + 
    8*Lambda1201*ms2)*(2*(Lambda1201*(Lambda1200^2 + Lambda1201^2) + 
      Lambda1211*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)) + 
    Lambda1201*((-3*g1^2)/5 - 3*g2^2 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + 
        Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + 
        Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
        Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QH1p^2 - 2*gN^2*QH2p^2 - 
      2*gN^2*QSp^2)) + 
  (4*(2*Lambda1210*mH1I200 + Lambda1211*mH1I201 + Lambda1211*mH1I210) + 
    4*(Lambda1200*mH2I201 + Lambda1200*mH2I210 + 2*Lambda1210*mH2I211) + 
    8*Lambda1210*ms2)*
   (2*(Lambda1200*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + 
      Lambda1210*(Lambda1210^2 + Lambda1211^2)) + 
    Lambda1210*((-3*g1^2)/5 - 3*g2^2 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + 
        Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + 
        Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
        Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QH1p^2 - 2*gN^2*QH2p^2 - 
      2*gN^2*QSp^2)) + 
  (4*(Lambda1210*mH1I201 + Lambda1210*mH1I210 + 2*Lambda1211*mH1I211) + 
    4*(Lambda1201*mH2I201 + Lambda1201*mH2I210 + 2*Lambda1211*mH2I211) + 
    8*Lambda1211*ms2)*
   (2*(Lambda1201*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211) + 
      Lambda1211*(Lambda1210^2 + Lambda1211^2)) + 
    Lambda1211*((-3*g1^2)/5 - 3*g2^2 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + 
        Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + 
        Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
        Lambda1211^2) + 2*Lambdax^2 - 2*gN^2*QH1p^2 - 2*gN^2*QH2p^2 - 
      2*gN^2*QSp^2)) - 32*gN^4*MassBp^2*QSp^2*(9*Qdp^2 + 9*QDxbarp^2 + 
    9*QDxp^2 + 3*Qep^2 + 6*QH1p^2 + 6*QH2p^2 + 2*QHpbarp^2 + 2*QHpp^2 + 
    6*QLp^2 + 18*QQp^2 + 3*QSp^2 + 9*Qup^2) + 
  8*gN^3*QSp^3*(-4*gN*MassBp^2*QSp + gN*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)) + 
  4*gN^2*QHpbarp*QSp*((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 + 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QHpbarp^2 + 2*gN^2*QHpbarp*
     (3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
       QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + 
      (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 
      2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 
      2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 
      6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 
      3*(mu200 + mu211 + mu222)*Qup)) + 4*gN^2*QHpp*QSp*
   ((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 - 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QHpp^2 + 2*gN^2*QHpp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)) + 
  gN^3*(9*Qdp^2 + 9*QDxbarp^2 + 9*QDxp^2 + 3*Qep^2 + 6*QH1p^2 + 6*QH2p^2 + 
    2*QHpbarp^2 + 2*QHpp^2 + 6*QLp^2 + 18*QQp^2 + 3*QSp^2 + 9*Qup^2)*
   (-16*gN*MassBp^2*QSp^2 + 4*gN*QSp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)) + 
  (6*(Kappa00^2 + Kappa01^2 + Kappa02^2) + 6*gN^2*QDxp*QSp)*
   ((-8*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 + 
    (Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
    Kappa00*(Kappa00*mDx200 + Kappa10*mDx201 + Kappa20*mDx202) + 
    Kappa01*(Kappa01*mDx200 + Kappa11*mDx201 + Kappa21*mDx202) + 
    Kappa02*(Kappa02*mDx200 + Kappa12*mDx201 + Kappa22*mDx202) + 
    (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
    (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220 + 
    2*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + Kappa02*mDxbar220) + 
      Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
      Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + Kappa02*mDxbar222)) + 
    2*(Kappa00^2 + Kappa01^2 + Kappa02^2)*ms2 - 
    (2*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QDxp^2 + 2*gN^2*QDxp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TKappa00^2 + TKappa01^2 + TKappa02^2)) + 
  6*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
   ((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx201 + 
    Kappa10*(Kappa00*mDx200 + Kappa10*mDx201 + Kappa20*mDx202) + 
    Kappa11*(Kappa01*mDx200 + Kappa11*mDx201 + Kappa21*mDx202) + 
    Kappa12*(Kappa02*mDx200 + Kappa12*mDx201 + Kappa22*mDx202) + 
    (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx211 + 
    (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx221 + 
    2*(Kappa10*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + Kappa02*mDxbar220) + 
      Kappa11*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
      Kappa12*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + Kappa02*mDxbar222)) + 
    2*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*ms2 + 
    2*(TKappa00*TKappa10 + TKappa01*TKappa11 + TKappa02*TKappa12)) + 
  6*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*
   ((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx200 + 
    (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx210 + 
    Kappa00*(Kappa00*mDx210 + Kappa10*mDx211 + Kappa20*mDx212) + 
    Kappa01*(Kappa01*mDx210 + Kappa11*mDx211 + Kappa21*mDx212) + 
    Kappa02*(Kappa02*mDx210 + Kappa12*mDx211 + Kappa22*mDx212) + 
    (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx220 + 
    2*(Kappa00*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + Kappa12*mDxbar220) + 
      Kappa01*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
      Kappa02*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + Kappa12*mDxbar222)) + 
    2*(Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*ms2 + 
    2*(TKappa00*TKappa10 + TKappa01*TKappa11 + TKappa02*TKappa12)) + 
  (6*(Kappa10^2 + Kappa11^2 + Kappa12^2) + 6*gN^2*QDxp*QSp)*
   ((-8*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 + 
    (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 + 
    (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
    Kappa10*(Kappa00*mDx210 + Kappa10*mDx211 + Kappa20*mDx212) + 
    Kappa11*(Kappa01*mDx210 + Kappa11*mDx211 + Kappa21*mDx212) + 
    Kappa12*(Kappa02*mDx210 + Kappa12*mDx211 + Kappa22*mDx212) + 
    (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 + 
    2*(Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + Kappa12*mDxbar220) + 
      Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
      Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + Kappa12*mDxbar222)) + 
    2*(Kappa10^2 + Kappa11^2 + Kappa12^2)*ms2 - 
    (2*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QDxp^2 + 2*gN^2*QDxp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TKappa10^2 + TKappa11^2 + TKappa12^2)) + 
  (6*(Kappa00^2 + Kappa10^2 + Kappa20^2) + 6*gN^2*QDxbarp*QSp)*
   ((-8*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 + 
    2*(Kappa00*(Kappa00*mDx200 + Kappa10*mDx210 + Kappa20*mDx220) + 
      Kappa10*(Kappa00*mDx201 + Kappa10*mDx211 + Kappa20*mDx221) + 
      Kappa20*(Kappa00*mDx202 + Kappa10*mDx212 + Kappa20*mDx222)) + 
    (Kappa00^2 + Kappa10^2 + Kappa20^2)*mDxbar200 + 
    Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar201 + Kappa02*mDxbar202) + 
    Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar201 + Kappa12*mDxbar202) + 
    Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar201 + Kappa22*mDxbar202) + 
    (Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*mDxbar210 + 
    (Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*mDxbar220 + 
    2*(Kappa00^2 + Kappa10^2 + Kappa20^2)*ms2 + 
    (2*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QDxbarp^2 + 2*gN^2*QDxbarp*
     (3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
       QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + 
      (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 
      2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 
      2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 
      6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 
      3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TKappa00^2 + TKappa10^2 + TKappa20^2)) + 
  6*(Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*
   (2*(Kappa00*(Kappa01*mDx200 + Kappa11*mDx210 + Kappa21*mDx220) + 
      Kappa10*(Kappa01*mDx201 + Kappa11*mDx211 + Kappa21*mDx221) + 
      Kappa20*(Kappa01*mDx202 + Kappa11*mDx212 + Kappa21*mDx222)) + 
    (Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*mDxbar200 + 
    (Kappa01^2 + Kappa11^2 + Kappa21^2)*mDxbar210 + 
    Kappa00*(Kappa00*mDxbar210 + Kappa01*mDxbar211 + Kappa02*mDxbar212) + 
    Kappa10*(Kappa10*mDxbar210 + Kappa11*mDxbar211 + Kappa12*mDxbar212) + 
    Kappa20*(Kappa20*mDxbar210 + Kappa21*mDxbar211 + Kappa22*mDxbar212) + 
    (Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*mDxbar220 + 
    2*(Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*ms2 + 
    2*(TKappa00*TKappa01 + TKappa10*TKappa11 + TKappa20*TKappa21)) + 
  6*(Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*
   (2*(Kappa01*(Kappa00*mDx200 + Kappa10*mDx210 + Kappa20*mDx220) + 
      Kappa11*(Kappa00*mDx201 + Kappa10*mDx211 + Kappa20*mDx221) + 
      Kappa21*(Kappa00*mDx202 + Kappa10*mDx212 + Kappa20*mDx222)) + 
    (Kappa00^2 + Kappa10^2 + Kappa20^2)*mDxbar201 + 
    Kappa01*(Kappa00*mDxbar200 + Kappa01*mDxbar201 + Kappa02*mDxbar202) + 
    Kappa11*(Kappa10*mDxbar200 + Kappa11*mDxbar201 + Kappa12*mDxbar202) + 
    Kappa21*(Kappa20*mDxbar200 + Kappa21*mDxbar201 + Kappa22*mDxbar202) + 
    (Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*mDxbar211 + 
    (Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*mDxbar221 + 
    2*(Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*ms2 + 
    2*(TKappa00*TKappa01 + TKappa10*TKappa11 + TKappa20*TKappa21)) + 
  (6*(Kappa01^2 + Kappa11^2 + Kappa21^2) + 6*gN^2*QDxbarp*QSp)*
   ((-8*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 + 
    2*(Kappa01*(Kappa01*mDx200 + Kappa11*mDx210 + Kappa21*mDx220) + 
      Kappa11*(Kappa01*mDx201 + Kappa11*mDx211 + Kappa21*mDx221) + 
      Kappa21*(Kappa01*mDx202 + Kappa11*mDx212 + Kappa21*mDx222)) + 
    (Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*mDxbar201 + 
    (Kappa01^2 + Kappa11^2 + Kappa21^2)*mDxbar211 + 
    Kappa01*(Kappa00*mDxbar210 + Kappa01*mDxbar211 + Kappa02*mDxbar212) + 
    Kappa11*(Kappa10*mDxbar210 + Kappa11*mDxbar211 + Kappa12*mDxbar212) + 
    Kappa21*(Kappa20*mDxbar210 + Kappa21*mDxbar211 + Kappa22*mDxbar212) + 
    (Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*mDxbar221 + 
    2*(Kappa01^2 + Kappa11^2 + Kappa21^2)*ms2 + 
    (2*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QDxbarp^2 + 2*gN^2*QDxbarp*
     (3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
       QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + 
      (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 
      2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 
      2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 
      6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 
      3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TKappa01^2 + TKappa11^2 + TKappa21^2)) + 
  6*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
   ((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx202 + 
    Kappa20*(Kappa00*mDx200 + Kappa10*mDx201 + Kappa20*mDx202) + 
    Kappa21*(Kappa01*mDx200 + Kappa11*mDx201 + Kappa21*mDx202) + 
    Kappa22*(Kappa02*mDx200 + Kappa12*mDx201 + Kappa22*mDx202) + 
    (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx212 + 
    (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx222 + 
    2*(Kappa20*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + Kappa02*mDxbar220) + 
      Kappa21*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
      Kappa22*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + Kappa02*mDxbar222)) + 
    2*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*ms2 + 
    2*(TKappa00*TKappa20 + TKappa01*TKappa21 + TKappa02*TKappa22)) + 
  6*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*
   ((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx200 + 
    (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx210 + 
    (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx220 + 
    Kappa00*(Kappa00*mDx220 + Kappa10*mDx221 + Kappa20*mDx222) + 
    Kappa01*(Kappa01*mDx220 + Kappa11*mDx221 + Kappa21*mDx222) + 
    Kappa02*(Kappa02*mDx220 + Kappa12*mDx221 + Kappa22*mDx222) + 
    2*(Kappa00*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + Kappa22*mDxbar220) + 
      Kappa01*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
      Kappa02*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + Kappa22*mDxbar222)) + 
    2*(Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*ms2 + 
    2*(TKappa00*TKappa20 + TKappa01*TKappa21 + TKappa02*TKappa22)) + 
  6*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
   ((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx202 + 
    (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx212 + 
    Kappa20*(Kappa00*mDx210 + Kappa10*mDx211 + Kappa20*mDx212) + 
    Kappa21*(Kappa01*mDx210 + Kappa11*mDx211 + Kappa21*mDx212) + 
    Kappa22*(Kappa02*mDx210 + Kappa12*mDx211 + Kappa22*mDx212) + 
    (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx222 + 
    2*(Kappa20*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + Kappa12*mDxbar220) + 
      Kappa21*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
      Kappa22*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + Kappa12*mDxbar222)) + 
    2*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*ms2 + 
    2*(TKappa10*TKappa20 + TKappa11*TKappa21 + TKappa12*TKappa22)) + 
  6*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*
   ((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx201 + 
    (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx211 + 
    (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx221 + 
    Kappa10*(Kappa00*mDx220 + Kappa10*mDx221 + Kappa20*mDx222) + 
    Kappa11*(Kappa01*mDx220 + Kappa11*mDx221 + Kappa21*mDx222) + 
    Kappa12*(Kappa02*mDx220 + Kappa12*mDx221 + Kappa22*mDx222) + 
    2*(Kappa10*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + Kappa22*mDxbar220) + 
      Kappa11*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
      Kappa12*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + Kappa22*mDxbar222)) + 
    2*(Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*ms2 + 
    2*(TKappa10*TKappa20 + TKappa11*TKappa21 + TKappa12*TKappa22)) + 
  6*(Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*
   (2*(Kappa02*(Kappa00*mDx200 + Kappa10*mDx210 + Kappa20*mDx220) + 
      Kappa12*(Kappa00*mDx201 + Kappa10*mDx211 + Kappa20*mDx221) + 
      Kappa22*(Kappa00*mDx202 + Kappa10*mDx212 + Kappa20*mDx222)) + 
    (Kappa00^2 + Kappa10^2 + Kappa20^2)*mDxbar202 + 
    Kappa02*(Kappa00*mDxbar200 + Kappa01*mDxbar201 + Kappa02*mDxbar202) + 
    Kappa12*(Kappa10*mDxbar200 + Kappa11*mDxbar201 + Kappa12*mDxbar202) + 
    Kappa22*(Kappa20*mDxbar200 + Kappa21*mDxbar201 + Kappa22*mDxbar202) + 
    (Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*mDxbar212 + 
    (Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*mDxbar222 + 
    2*(Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*ms2 + 
    2*(TKappa00*TKappa02 + TKappa10*TKappa12 + TKappa20*TKappa22)) + 
  6*(Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*
   (2*(Kappa00*(Kappa02*mDx200 + Kappa12*mDx210 + Kappa22*mDx220) + 
      Kappa10*(Kappa02*mDx201 + Kappa12*mDx211 + Kappa22*mDx221) + 
      Kappa20*(Kappa02*mDx202 + Kappa12*mDx212 + Kappa22*mDx222)) + 
    (Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*mDxbar200 + 
    (Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*mDxbar210 + 
    (Kappa02^2 + Kappa12^2 + Kappa22^2)*mDxbar220 + 
    Kappa00*(Kappa00*mDxbar220 + Kappa01*mDxbar221 + Kappa02*mDxbar222) + 
    Kappa10*(Kappa10*mDxbar220 + Kappa11*mDxbar221 + Kappa12*mDxbar222) + 
    Kappa20*(Kappa20*mDxbar220 + Kappa21*mDxbar221 + Kappa22*mDxbar222) + 
    2*(Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*ms2 + 
    2*(TKappa00*TKappa02 + TKappa10*TKappa12 + TKappa20*TKappa22)) + 
  6*(Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*
   (2*(Kappa02*(Kappa01*mDx200 + Kappa11*mDx210 + Kappa21*mDx220) + 
      Kappa12*(Kappa01*mDx201 + Kappa11*mDx211 + Kappa21*mDx221) + 
      Kappa22*(Kappa01*mDx202 + Kappa11*mDx212 + Kappa21*mDx222)) + 
    (Kappa00*Kappa01 + Kappa10*Kappa11 + Kappa20*Kappa21)*mDxbar202 + 
    (Kappa01^2 + Kappa11^2 + Kappa21^2)*mDxbar212 + 
    Kappa02*(Kappa00*mDxbar210 + Kappa01*mDxbar211 + Kappa02*mDxbar212) + 
    Kappa12*(Kappa10*mDxbar210 + Kappa11*mDxbar211 + Kappa12*mDxbar212) + 
    Kappa22*(Kappa20*mDxbar210 + Kappa21*mDxbar211 + Kappa22*mDxbar212) + 
    (Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*mDxbar222 + 
    2*(Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*ms2 + 
    2*(TKappa01*TKappa02 + TKappa11*TKappa12 + TKappa21*TKappa22)) + 
  6*(Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*
   (2*(Kappa01*(Kappa02*mDx200 + Kappa12*mDx210 + Kappa22*mDx220) + 
      Kappa11*(Kappa02*mDx201 + Kappa12*mDx211 + Kappa22*mDx221) + 
      Kappa21*(Kappa02*mDx202 + Kappa12*mDx212 + Kappa22*mDx222)) + 
    (Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*mDxbar201 + 
    (Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*mDxbar211 + 
    (Kappa02^2 + Kappa12^2 + Kappa22^2)*mDxbar221 + 
    Kappa01*(Kappa00*mDxbar220 + Kappa01*mDxbar221 + Kappa02*mDxbar222) + 
    Kappa11*(Kappa10*mDxbar220 + Kappa11*mDxbar221 + Kappa12*mDxbar222) + 
    Kappa21*(Kappa20*mDxbar220 + Kappa21*mDxbar221 + Kappa22*mDxbar222) + 
    2*(Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*ms2 + 
    2*(TKappa01*TKappa02 + TKappa11*TKappa12 + TKappa21*TKappa22)) + 
  (6*(Kappa02^2 + Kappa12^2 + Kappa22^2) + 6*gN^2*QDxbarp*QSp)*
   ((-8*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 + 
    2*(Kappa02*(Kappa02*mDx200 + Kappa12*mDx210 + Kappa22*mDx220) + 
      Kappa12*(Kappa02*mDx201 + Kappa12*mDx211 + Kappa22*mDx221) + 
      Kappa22*(Kappa02*mDx202 + Kappa12*mDx212 + Kappa22*mDx222)) + 
    (Kappa00*Kappa02 + Kappa10*Kappa12 + Kappa20*Kappa22)*mDxbar202 + 
    (Kappa01*Kappa02 + Kappa11*Kappa12 + Kappa21*Kappa22)*mDxbar212 + 
    (Kappa02^2 + Kappa12^2 + Kappa22^2)*mDxbar222 + 
    Kappa02*(Kappa00*mDxbar220 + Kappa01*mDxbar221 + Kappa02*mDxbar222) + 
    Kappa12*(Kappa10*mDxbar220 + Kappa11*mDxbar221 + Kappa12*mDxbar222) + 
    Kappa22*(Kappa20*mDxbar220 + Kappa21*mDxbar221 + Kappa22*mDxbar222) + 
    2*(Kappa02^2 + Kappa12^2 + Kappa22^2)*ms2 + 
    (2*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QDxbarp^2 + 2*gN^2*QDxbarp*
     (3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
       QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + 
      (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 
      2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 
      2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 
      6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 
      3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TKappa02^2 + TKappa12^2 + TKappa22^2)) + 
  (6*(Kappa20^2 + Kappa21^2 + Kappa22^2) + 6*gN^2*QDxp*QSp)*
   ((-8*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 + 
    (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx202 + 
    (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx212 + 
    (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222 + 
    Kappa20*(Kappa00*mDx220 + Kappa10*mDx221 + Kappa20*mDx222) + 
    Kappa21*(Kappa01*mDx220 + Kappa11*mDx221 + Kappa21*mDx222) + 
    Kappa22*(Kappa02*mDx220 + Kappa12*mDx221 + Kappa22*mDx222) + 
    2*(Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + Kappa22*mDxbar220) + 
      Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
      Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + Kappa22*mDxbar222)) + 
    2*(Kappa20^2 + Kappa21^2 + Kappa22^2)*ms2 - 
    (2*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QDxp^2 + 2*gN^2*QDxp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TKappa20^2 + TKappa21^2 + TKappa22^2)) + 
  (4*(Lambda1200^2 + Lambda1201^2) + 4*gN^2*QH2p*QSp)*
   ((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 + 
    2*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I210) + 
      Lambda1201*(Lambda1200*mH1I201 + Lambda1201*mH1I211)) + 
    (Lambda1200^2 + Lambda1201^2)*mH2I200 + 
    Lambda1200*(Lambda1200*mH2I200 + Lambda1210*mH2I201) + 
    Lambda1201*(Lambda1201*mH2I200 + Lambda1211*mH2I201) + 
    (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + 
    2*(Lambda1200^2 + Lambda1201^2)*ms2 + 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QH2p^2 + 2*gN^2*QH2p*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TLambda1200^2 + TLambda1201^2)) + 
  (4*(Lambda1200^2 + Lambda1210^2) + 4*gN^2*QH1p*QSp)*
   ((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 + (Lambda1200^2 + Lambda1210^2)*
     mH1I200 + Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
    Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + 
    (Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*mH1I210 + 
    2*(Lambda1200*(Lambda1200*mH2I200 + Lambda1210*mH2I210) + 
      Lambda1210*(Lambda1200*mH2I201 + Lambda1210*mH2I211)) + 
    2*(Lambda1200^2 + Lambda1210^2)*ms2 - 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QH1p^2 + 2*gN^2*QH1p*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TLambda1200^2 + TLambda1210^2)) + 
  4*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*
   (2*(Lambda1210*(Lambda1200*mH1I200 + Lambda1201*mH1I210) + 
      Lambda1211*(Lambda1200*mH1I201 + Lambda1201*mH1I211)) + 
    (Lambda1200^2 + Lambda1201^2)*mH2I201 + 
    Lambda1210*(Lambda1200*mH2I200 + Lambda1210*mH2I201) + 
    Lambda1211*(Lambda1201*mH2I200 + Lambda1211*mH2I201) + 
    (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I211 + 
    2*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*ms2 + 
    2*(TLambda1200*TLambda1210 + TLambda1201*TLambda1211)) + 
  4*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*
   (2*(Lambda1200*(Lambda1210*mH1I200 + Lambda1211*mH1I210) + 
      Lambda1201*(Lambda1210*mH1I201 + Lambda1211*mH1I211)) + 
    (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I200 + 
    (Lambda1210^2 + Lambda1211^2)*mH2I210 + 
    Lambda1200*(Lambda1200*mH2I210 + Lambda1210*mH2I211) + 
    Lambda1201*(Lambda1201*mH2I210 + Lambda1211*mH2I211) + 
    2*(Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*ms2 + 
    2*(TLambda1200*TLambda1210 + TLambda1201*TLambda1211)) + 
  4*(Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*
   ((Lambda1200^2 + Lambda1210^2)*mH1I201 + 
    Lambda1201*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
    Lambda1211*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + 
    (Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*mH1I211 + 
    2*(Lambda1201*(Lambda1200*mH2I200 + Lambda1210*mH2I210) + 
      Lambda1211*(Lambda1200*mH2I201 + Lambda1210*mH2I211)) + 
    2*(Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*ms2 + 
    2*(TLambda1200*TLambda1201 + TLambda1210*TLambda1211)) + 
  4*(Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*
   ((Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*mH1I200 + 
    (Lambda1201^2 + Lambda1211^2)*mH1I210 + 
    Lambda1200*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
    Lambda1210*(Lambda1210*mH1I210 + Lambda1211*mH1I211) + 
    2*(Lambda1200*(Lambda1201*mH2I200 + Lambda1211*mH2I210) + 
      Lambda1210*(Lambda1201*mH2I201 + Lambda1211*mH2I211)) + 
    2*(Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*ms2 + 
    2*(TLambda1200*TLambda1201 + TLambda1210*TLambda1211)) + 
  (4*(Lambda1201^2 + Lambda1211^2) + 4*gN^2*QH1p*QSp)*
   ((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 + 
    (Lambda1200*Lambda1201 + Lambda1210*Lambda1211)*mH1I201 + 
    (Lambda1201^2 + Lambda1211^2)*mH1I211 + 
    Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
    Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211) + 
    2*(Lambda1201*(Lambda1201*mH2I200 + Lambda1211*mH2I210) + 
      Lambda1211*(Lambda1201*mH2I201 + Lambda1211*mH2I211)) + 
    2*(Lambda1201^2 + Lambda1211^2)*ms2 - 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QH1p^2 + 2*gN^2*QH1p*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TLambda1201^2 + TLambda1211^2)) + 
  (4*(Lambda1210^2 + Lambda1211^2) + 4*gN^2*QH2p*QSp)*
   ((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 + 
    2*(Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I210) + 
      Lambda1211*(Lambda1210*mH1I201 + Lambda1211*mH1I211)) + 
    (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
    (Lambda1210^2 + Lambda1211^2)*mH2I211 + 
    Lambda1210*(Lambda1200*mH2I210 + Lambda1210*mH2I211) + 
    Lambda1211*(Lambda1201*mH2I210 + Lambda1211*mH2I211) + 
    2*(Lambda1210^2 + Lambda1211^2)*ms2 + 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QH2p^2 + 2*gN^2*QH2p*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TLambda1210^2 + TLambda1211^2)) + 
  (6*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + 
      Kappa20^2 + Kappa21^2 + Kappa22^2) + 
    4*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2) + 
    4*Lambdax^2 + 2*gN^2*QSp^2)*
   (6*((Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx202 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
      (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx212 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 + 
      (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222) + 
    6*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + Kappa02*mDxbar220) + 
      Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + Kappa12*mDxbar220) + 
      Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + Kappa22*mDxbar220) + 
      Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + Kappa02*mDxbar221) + 
      Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + Kappa12*mDxbar221) + 
      Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + Kappa22*mDxbar221) + 
      Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + Kappa02*mDxbar222) + 
      Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + Kappa12*mDxbar222) + 
      Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + Kappa22*mDxbar222)) + 
    4*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
      Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + 
      Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
      Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211)) + 
    4*((Lambda1200^2 + Lambda1201^2)*mH2I200 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + 
      (Lambda1210^2 + Lambda1211^2)*mH2I211) + 
    6*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
      Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*ms2 + 
    4*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*ms2 + 
    4*Lambdax^2*(mHd2 + mHu2 + ms2) - 8*gN^2*MassBp^2*QSp^2 + 
    2*gN^2*QSp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    6*(TKappa00^2 + TKappa01^2 + TKappa02^2 + TKappa10^2 + TKappa11^2 + 
      TKappa12^2 + TKappa20^2 + TKappa21^2 + TKappa22^2) + 
    4*(TLambda1200^2 + TLambda1201^2 + TLambda1210^2 + TLambda1211^2) + 
    4*TLambdax^2) + 12*TKappa00*((-4*g1^2*TKappa00)/15 - 
    (16*g3^2*TKappa00)/3 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa00 + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*TKappa00 + 
    2*Lambdax^2*TKappa00 - 2*gN^2*QDxbarp^2*TKappa00 - 
    2*gN^2*QDxp^2*TKappa00 - 2*gN^2*QSp^2*TKappa00 + 
    3*(Kappa00*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02) + 
      Kappa10*(Kappa10*TKappa00 + Kappa11*TKappa01 + Kappa12*TKappa02) + 
      Kappa20*(Kappa20*TKappa00 + Kappa21*TKappa01 + Kappa22*TKappa02)) + 
    3*((Kappa00^2 + Kappa01^2 + Kappa02^2)*TKappa00 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa10 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa20) + 
    Kappa00*((8*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      4*gN^2*MassBp*QDxbarp^2 + 4*gN^2*MassBp*QDxp^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 12*TKappa01*((-4*g1^2*TKappa01)/15 - 
    (16*g3^2*TKappa01)/3 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa01 + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*TKappa01 + 
    2*Lambdax^2*TKappa01 - 2*gN^2*QDxbarp^2*TKappa01 - 
    2*gN^2*QDxp^2*TKappa01 - 2*gN^2*QSp^2*TKappa01 + 
    3*(Kappa01*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02) + 
      Kappa11*(Kappa10*TKappa00 + Kappa11*TKappa01 + Kappa12*TKappa02) + 
      Kappa21*(Kappa20*TKappa00 + Kappa21*TKappa01 + Kappa22*TKappa02)) + 
    3*((Kappa00^2 + Kappa01^2 + Kappa02^2)*TKappa01 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa11 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa21) + 
    Kappa01*((8*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      4*gN^2*MassBp*QDxbarp^2 + 4*gN^2*MassBp*QDxp^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 12*TKappa02*((-4*g1^2*TKappa02)/15 - 
    (16*g3^2*TKappa02)/3 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa02 + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*TKappa02 + 
    2*Lambdax^2*TKappa02 - 2*gN^2*QDxbarp^2*TKappa02 - 
    2*gN^2*QDxp^2*TKappa02 - 2*gN^2*QSp^2*TKappa02 + 
    3*(Kappa02*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02) + 
      Kappa12*(Kappa10*TKappa00 + Kappa11*TKappa01 + Kappa12*TKappa02) + 
      Kappa22*(Kappa20*TKappa00 + Kappa21*TKappa01 + Kappa22*TKappa02)) + 
    3*((Kappa00^2 + Kappa01^2 + Kappa02^2)*TKappa02 + 
      (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa12 + 
      (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa22) + 
    Kappa02*((8*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      4*gN^2*MassBp*QDxbarp^2 + 4*gN^2*MassBp*QDxp^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 12*TKappa10*((-4*g1^2*TKappa10)/15 - 
    (16*g3^2*TKappa10)/3 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa10 + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*TKappa10 + 
    2*Lambdax^2*TKappa10 - 2*gN^2*QDxbarp^2*TKappa10 - 
    2*gN^2*QDxp^2*TKappa10 - 2*gN^2*QSp^2*TKappa10 + 
    3*(Kappa00*(Kappa00*TKappa10 + Kappa01*TKappa11 + Kappa02*TKappa12) + 
      Kappa10*(Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12) + 
      Kappa20*(Kappa20*TKappa10 + Kappa21*TKappa11 + Kappa22*TKappa12)) + 
    3*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa00 + 
      (Kappa10^2 + Kappa11^2 + Kappa12^2)*TKappa10 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa20) + 
    Kappa10*((8*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      4*gN^2*MassBp*QDxbarp^2 + 4*gN^2*MassBp*QDxp^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 12*TKappa11*((-4*g1^2*TKappa11)/15 - 
    (16*g3^2*TKappa11)/3 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa11 + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*TKappa11 + 
    2*Lambdax^2*TKappa11 - 2*gN^2*QDxbarp^2*TKappa11 - 
    2*gN^2*QDxp^2*TKappa11 - 2*gN^2*QSp^2*TKappa11 + 
    3*(Kappa01*(Kappa00*TKappa10 + Kappa01*TKappa11 + Kappa02*TKappa12) + 
      Kappa11*(Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12) + 
      Kappa21*(Kappa20*TKappa10 + Kappa21*TKappa11 + Kappa22*TKappa12)) + 
    3*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa01 + 
      (Kappa10^2 + Kappa11^2 + Kappa12^2)*TKappa11 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa21) + 
    Kappa11*((8*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      4*gN^2*MassBp*QDxbarp^2 + 4*gN^2*MassBp*QDxp^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 12*TKappa12*((-4*g1^2*TKappa12)/15 - 
    (16*g3^2*TKappa12)/3 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa12 + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*TKappa12 + 
    2*Lambdax^2*TKappa12 - 2*gN^2*QDxbarp^2*TKappa12 - 
    2*gN^2*QDxp^2*TKappa12 - 2*gN^2*QSp^2*TKappa12 + 
    3*(Kappa02*(Kappa00*TKappa10 + Kappa01*TKappa11 + Kappa02*TKappa12) + 
      Kappa12*(Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12) + 
      Kappa22*(Kappa20*TKappa10 + Kappa21*TKappa11 + Kappa22*TKappa12)) + 
    3*((Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*TKappa02 + 
      (Kappa10^2 + Kappa11^2 + Kappa12^2)*TKappa12 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa22) + 
    Kappa12*((8*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      4*gN^2*MassBp*QDxbarp^2 + 4*gN^2*MassBp*QDxp^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 12*TKappa20*((-4*g1^2*TKappa20)/15 - 
    (16*g3^2*TKappa20)/3 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa20 + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*TKappa20 + 
    2*Lambdax^2*TKappa20 - 2*gN^2*QDxbarp^2*TKappa20 - 
    2*gN^2*QDxp^2*TKappa20 - 2*gN^2*QSp^2*TKappa20 + 
    3*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa00 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa10 + 
      (Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa20) + 
    3*(Kappa00*(Kappa00*TKappa20 + Kappa01*TKappa21 + Kappa02*TKappa22) + 
      Kappa10*(Kappa10*TKappa20 + Kappa11*TKappa21 + Kappa12*TKappa22) + 
      Kappa20*(Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)) + 
    Kappa20*((8*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      4*gN^2*MassBp*QDxbarp^2 + 4*gN^2*MassBp*QDxp^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 12*TKappa21*((-4*g1^2*TKappa21)/15 - 
    (16*g3^2*TKappa21)/3 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa21 + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*TKappa21 + 
    2*Lambdax^2*TKappa21 - 2*gN^2*QDxbarp^2*TKappa21 - 
    2*gN^2*QDxp^2*TKappa21 - 2*gN^2*QSp^2*TKappa21 + 
    3*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa01 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa11 + 
      (Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa21) + 
    3*(Kappa01*(Kappa00*TKappa20 + Kappa01*TKappa21 + Kappa02*TKappa22) + 
      Kappa11*(Kappa10*TKappa20 + Kappa11*TKappa21 + Kappa12*TKappa22) + 
      Kappa21*(Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)) + 
    Kappa21*((8*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      4*gN^2*MassBp*QDxbarp^2 + 4*gN^2*MassBp*QDxp^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 12*TKappa22*((-4*g1^2*TKappa22)/15 - 
    (16*g3^2*TKappa22)/3 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa22 + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*TKappa22 + 
    2*Lambdax^2*TKappa22 - 2*gN^2*QDxbarp^2*TKappa22 - 
    2*gN^2*QDxp^2*TKappa22 - 2*gN^2*QSp^2*TKappa22 + 
    3*((Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*TKappa02 + 
      (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*TKappa12 + 
      (Kappa20^2 + Kappa21^2 + Kappa22^2)*TKappa22) + 
    3*(Kappa02*(Kappa00*TKappa20 + Kappa01*TKappa21 + Kappa02*TKappa22) + 
      Kappa12*(Kappa10*TKappa20 + Kappa11*TKappa21 + Kappa12*TKappa22) + 
      Kappa22*(Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22)) + 
    Kappa22*((8*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      4*gN^2*MassBp*QDxbarp^2 + 4*gN^2*MassBp*QDxp^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 8*TLambda1200*((-3*g1^2*TLambda1200)/5 - 
    3*g2^2*TLambda1200 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*
     TLambda1200 + 2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
      Lambda1211^2)*TLambda1200 + 2*Lambdax^2*TLambda1200 - 
    2*gN^2*QH1p^2*TLambda1200 - 2*gN^2*QH2p^2*TLambda1200 - 
    2*gN^2*QSp^2*TLambda1200 + 
    3*(Lambda1200*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201) + 
      Lambda1210*(Lambda1210*TLambda1200 + Lambda1211*TLambda1201)) + 
    3*((Lambda1200^2 + Lambda1201^2)*TLambda1200 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1210) + 
    Lambda1200*((6*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*QH1p^2 + 
      4*gN^2*MassBp*QH2p^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 8*TLambda1201*((-3*g1^2*TLambda1201)/5 - 
    3*g2^2*TLambda1201 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*
     TLambda1201 + 2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
      Lambda1211^2)*TLambda1201 + 2*Lambdax^2*TLambda1201 - 
    2*gN^2*QH1p^2*TLambda1201 - 2*gN^2*QH2p^2*TLambda1201 - 
    2*gN^2*QSp^2*TLambda1201 + 
    3*(Lambda1201*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201) + 
      Lambda1211*(Lambda1210*TLambda1200 + Lambda1211*TLambda1201)) + 
    3*((Lambda1200^2 + Lambda1201^2)*TLambda1201 + 
      (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1211) + 
    Lambda1201*((6*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*QH1p^2 + 
      4*gN^2*MassBp*QH2p^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 8*TLambda1210*((-3*g1^2*TLambda1210)/5 - 
    3*g2^2*TLambda1210 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*
     TLambda1210 + 2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
      Lambda1211^2)*TLambda1210 + 2*Lambdax^2*TLambda1210 - 
    2*gN^2*QH1p^2*TLambda1210 - 2*gN^2*QH2p^2*TLambda1210 - 
    2*gN^2*QSp^2*TLambda1210 + 
    3*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1200 + 
      (Lambda1210^2 + Lambda1211^2)*TLambda1210) + 
    3*(Lambda1200*(Lambda1200*TLambda1210 + Lambda1201*TLambda1211) + 
      Lambda1210*(Lambda1210*TLambda1210 + Lambda1211*TLambda1211)) + 
    Lambda1210*((6*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*QH1p^2 + 
      4*gN^2*MassBp*QH2p^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 8*TLambda1211*((-3*g1^2*TLambda1211)/5 - 
    3*g2^2*TLambda1211 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + 
      Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*
     TLambda1211 + 2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
      Lambda1211^2)*TLambda1211 + 2*Lambdax^2*TLambda1211 - 
    2*gN^2*QH1p^2*TLambda1211 - 2*gN^2*QH2p^2*TLambda1211 - 
    2*gN^2*QSp^2*TLambda1211 + 
    3*((Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*TLambda1201 + 
      (Lambda1210^2 + Lambda1211^2)*TLambda1211) + 
    3*(Lambda1201*(Lambda1200*TLambda1210 + Lambda1201*TLambda1211) + 
      Lambda1211*(Lambda1210*TLambda1210 + Lambda1211*TLambda1211)) + 
    Lambda1211*((6*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*QH1p^2 + 
      4*gN^2*MassBp*QH2p^2 + 4*gN^2*MassBp*QSp^2 + 
      6*(Kappa00*TKappa00 + Kappa01*TKappa01 + Kappa02*TKappa02 + 
        Kappa10*TKappa10 + Kappa11*TKappa11 + Kappa12*TKappa12 + 
        Kappa20*TKappa20 + Kappa21*TKappa21 + Kappa22*TKappa22) + 
      4*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
        Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
      4*Lambdax*TLambdax)) + 6*gN^2*Qdp*QSp*((-8*g1^2*MassB^2)/15 - 
    (32*g3^2*MassG^2)/3 + (2*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - 
       mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - 
       mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - 
       ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 
       2*(mu200 + mu211 + mu222)))/5 - 8*gN^2*MassBp^2*Qdp^2 + 
    2*gN^2*Qdp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    4*(TYd00^2 + TYd01^2 + TYd02^2) + 4*mHd2*(Yd00^2 + Yd01^2 + Yd02^2) + 
    4*(Yd00*(mq200*Yd00 + mq210*Yd01 + mq220*Yd02) + 
      Yd01*(mq201*Yd00 + mq211*Yd01 + mq221*Yd02) + 
      Yd02*(mq202*Yd00 + mq212*Yd01 + mq222*Yd02)) + 
    2*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
      Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
      Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22)) + 
    2*(md200*(Yd00^2 + Yd01^2 + Yd02^2) + md210*(Yd00*Yd10 + Yd01*Yd11 + 
        Yd02*Yd12) + md220*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22))) + 
  6*gN^2*Qdp*QSp*((-8*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 + 
    (2*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*Qdp^2 + 2*gN^2*Qdp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    4*(TYd10^2 + TYd11^2 + TYd12^2) + 4*mHd2*(Yd10^2 + Yd11^2 + Yd12^2) + 
    4*(Yd10*(mq200*Yd10 + mq210*Yd11 + mq220*Yd12) + 
      Yd11*(mq201*Yd10 + mq211*Yd11 + mq221*Yd12) + 
      Yd12*(mq202*Yd10 + mq212*Yd11 + mq222*Yd12)) + 
    2*(Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
      Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
      Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22)) + 
    2*(md201*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      md211*(Yd10^2 + Yd11^2 + Yd12^2) + md221*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22))) + 6*gN^2*Qdp*QSp*((-8*g1^2*MassB^2)/15 - 
    (32*g3^2*MassG^2)/3 + (2*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - 
       mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - 
       mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - 
       ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 
       2*(mu200 + mu211 + mu222)))/5 - 8*gN^2*MassBp^2*Qdp^2 + 
    2*gN^2*Qdp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    4*(TYd20^2 + TYd21^2 + TYd22^2) + 4*mHd2*(Yd20^2 + Yd21^2 + Yd22^2) + 
    2*(Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
      Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
      Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
    4*(Yd20*(mq200*Yd20 + mq210*Yd21 + mq220*Yd22) + 
      Yd21*(mq201*Yd20 + mq211*Yd21 + mq221*Yd22) + 
      Yd22*(mq202*Yd20 + mq212*Yd21 + mq222*Yd22)) + 
    2*(md202*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      md212*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      md222*(Yd20^2 + Yd21^2 + Yd22^2))) + 
  4*gN^2*QLp*QSp*((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 - 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QLp^2 + 2*gN^2*QLp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TYe00^2 + TYe10^2 + TYe20^2) + Ye00*(ml200*Ye00 + ml201*Ye01 + 
      ml202*Ye02) + Ye10*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
    2*mHd2*(Ye00^2 + Ye10^2 + Ye20^2) + ml200*(Ye00^2 + Ye10^2 + Ye20^2) + 
    2*(Ye00*(me200*Ye00 + me210*Ye10 + me220*Ye20) + 
      Ye10*(me201*Ye00 + me211*Ye10 + me221*Ye20) + 
      Ye20*(me202*Ye00 + me212*Ye10 + me222*Ye20)) + 
    ml210*(Ye00*Ye01 + Ye10*Ye11 + Ye20*Ye21) + 
    Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22) + 
    ml220*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22)) + 
  4*gN^2*QLp*QSp*((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 - 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QLp^2 + 2*gN^2*QLp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TYe01^2 + TYe11^2 + TYe21^2) + Ye01*(ml210*Ye00 + ml211*Ye01 + 
      ml212*Ye02) + Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
    ml201*(Ye00*Ye01 + Ye10*Ye11 + Ye20*Ye21) + 
    2*mHd2*(Ye01^2 + Ye11^2 + Ye21^2) + ml211*(Ye01^2 + Ye11^2 + Ye21^2) + 
    2*(Ye01*(me200*Ye01 + me210*Ye11 + me220*Ye21) + 
      Ye11*(me201*Ye01 + me211*Ye11 + me221*Ye21) + 
      Ye21*(me202*Ye01 + me212*Ye11 + me222*Ye21)) + 
    Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + 
    ml221*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22)) + 
  4*gN^2*QLp*QSp*((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 - 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QLp^2 + 2*gN^2*QLp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TYe02^2 + TYe12^2 + TYe22^2) + Ye02*(ml220*Ye00 + ml221*Ye01 + 
      ml222*Ye02) + Ye12*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
    Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22) + 
    ml202*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22) + 
    ml212*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22) + 
    2*mHd2*(Ye02^2 + Ye12^2 + Ye22^2) + ml222*(Ye02^2 + Ye12^2 + Ye22^2) + 
    2*(Ye02*(me200*Ye02 + me210*Ye12 + me220*Ye22) + 
      Ye12*(me201*Ye02 + me211*Ye12 + me221*Ye22) + 
      Ye22*(me202*Ye02 + me212*Ye12 + me222*Ye22))) + 
  (4*Lambdax^2 + 4*gN^2*QH1p*QSp)*((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 + 
    2*Lambdax^2*mHd2 + 2*Lambdax^2*mHu2 + 2*Lambdax^2*ms2 - 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QH1p^2 + 2*gN^2*QH1p*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*TLambdax^2 + 6*(TYd00^2 + TYd01^2 + TYd02^2 + TYd10^2 + TYd11^2 + 
      TYd12^2 + TYd20^2 + TYd21^2 + TYd22^2) + 
    2*(TYe00^2 + TYe01^2 + TYe02^2 + TYe10^2 + TYe11^2 + TYe12^2 + TYe20^2 + 
      TYe21^2 + TYe22^2) + 6*mHd2*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    6*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
      Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
      Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
      Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
      Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
      Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
      Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22) + 
      Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22) + 
      Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
    6*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
      Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
      Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
      Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
      Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
      Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
      Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + 
      Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + 
      Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) + 
    2*mHd2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + 2*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
      Ye10*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
      Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
      Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
      Ye11*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
      Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
      Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22) + 
      Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22) + 
      Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
    2*(Ye00*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
      Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
      Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
      Ye10*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
      Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
      Ye12*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
      Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22) + 
      Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + 
      Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22))) + 
  2*gN^2*Qep*QSp*((-24*g1^2*MassB^2)/5 + 
    (6*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*Qep^2 + 2*gN^2*Qep*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    4*(TYe00^2 + TYe01^2 + TYe02^2) + 4*mHd2*(Ye00^2 + Ye01^2 + Ye02^2) + 
    4*(Ye00*(ml200*Ye00 + ml210*Ye01 + ml220*Ye02) + 
      Ye01*(ml201*Ye00 + ml211*Ye01 + ml221*Ye02) + 
      Ye02*(ml202*Ye00 + ml212*Ye01 + ml222*Ye02)) + 
    2*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
      Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
      Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22)) + 
    2*(me200*(Ye00^2 + Ye01^2 + Ye02^2) + me210*(Ye00*Ye10 + Ye01*Ye11 + 
        Ye02*Ye12) + me220*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22))) + 
  2*gN^2*Qep*QSp*((-24*g1^2*MassB^2)/5 + 
    (6*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*Qep^2 + 2*gN^2*Qep*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    4*(TYe10^2 + TYe11^2 + TYe12^2) + 4*mHd2*(Ye10^2 + Ye11^2 + Ye12^2) + 
    4*(Ye10*(ml200*Ye10 + ml210*Ye11 + ml220*Ye12) + 
      Ye11*(ml201*Ye10 + ml211*Ye11 + ml221*Ye12) + 
      Ye12*(ml202*Ye10 + ml212*Ye11 + ml222*Ye12)) + 
    2*(Ye10*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
      Ye11*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
      Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22)) + 
    2*(me201*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      me211*(Ye10^2 + Ye11^2 + Ye12^2) + me221*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22))) + 2*gN^2*Qep*QSp*((-24*g1^2*MassB^2)/5 + 
    (6*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*Qep^2 + 2*gN^2*Qep*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    4*(TYe20^2 + TYe21^2 + TYe22^2) + 4*mHd2*(Ye20^2 + Ye21^2 + Ye22^2) + 
    2*(Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
      Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
      Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
    4*(Ye20*(ml200*Ye20 + ml210*Ye21 + ml220*Ye22) + 
      Ye21*(ml201*Ye20 + ml211*Ye21 + ml221*Ye22) + 
      Ye22*(ml202*Ye20 + ml212*Ye21 + ml222*Ye22)) + 
    2*(me202*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      me212*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      me222*(Ye20^2 + Ye21^2 + Ye22^2))) + 
  12*gN^2*QQp*QSp*((-2*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 - 
    6*g2^2*MassWB^2 + (g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - 
       mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - 
       mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - 
       ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 
       2*(mu200 + mu211 + mu222)))/5 - 8*gN^2*MassBp^2*QQp^2 + 
    2*gN^2*QQp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TYd00^2 + TYd10^2 + TYd20^2) + 2*(TYu00^2 + TYu10^2 + TYu20^2) + 
    Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
    Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
    2*mHd2*(Yd00^2 + Yd10^2 + Yd20^2) + mq200*(Yd00^2 + Yd10^2 + Yd20^2) + 
    2*(Yd00*(md200*Yd00 + md210*Yd10 + md220*Yd20) + 
      Yd10*(md201*Yd00 + md211*Yd10 + md221*Yd20) + 
      Yd20*(md202*Yd00 + md212*Yd10 + md222*Yd20)) + 
    mq210*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21) + 
    Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + 
    mq220*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + 
    Yu00*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
    Yu10*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
    2*mHu2*(Yu00^2 + Yu10^2 + Yu20^2) + mq200*(Yu00^2 + Yu10^2 + Yu20^2) + 
    2*(Yu00*(mu200*Yu00 + mu210*Yu10 + mu220*Yu20) + 
      Yu10*(mu201*Yu00 + mu211*Yu10 + mu221*Yu20) + 
      Yu20*(mu202*Yu00 + mu212*Yu10 + mu222*Yu20)) + 
    mq210*(Yu00*Yu01 + Yu10*Yu11 + Yu20*Yu21) + 
    Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22) + 
    mq220*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22)) + 
  12*gN^2*QQp*QSp*((-2*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 - 
    6*g2^2*MassWB^2 + (g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - 
       mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - 
       mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - 
       ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 
       2*(mu200 + mu211 + mu222)))/5 - 8*gN^2*MassBp^2*QQp^2 + 
    2*gN^2*QQp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TYd01^2 + TYd11^2 + TYd21^2) + 2*(TYu01^2 + TYu11^2 + TYu21^2) + 
    Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
    Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
    mq201*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21) + 
    2*mHd2*(Yd01^2 + Yd11^2 + Yd21^2) + mq211*(Yd01^2 + Yd11^2 + Yd21^2) + 
    2*(Yd01*(md200*Yd01 + md210*Yd11 + md220*Yd21) + 
      Yd11*(md201*Yd01 + md211*Yd11 + md221*Yd21) + 
      Yd21*(md202*Yd01 + md212*Yd11 + md222*Yd21)) + 
    Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + 
    mq221*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + 
    Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
    Yu11*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
    mq201*(Yu00*Yu01 + Yu10*Yu11 + Yu20*Yu21) + 
    2*mHu2*(Yu01^2 + Yu11^2 + Yu21^2) + mq211*(Yu01^2 + Yu11^2 + Yu21^2) + 
    2*(Yu01*(mu200*Yu01 + mu210*Yu11 + mu220*Yu21) + 
      Yu11*(mu201*Yu01 + mu211*Yu11 + mu221*Yu21) + 
      Yu21*(mu202*Yu01 + mu212*Yu11 + mu222*Yu21)) + 
    Yu21*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22) + 
    mq221*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22)) + 
  8*Lambdax*(mHd2 + mHu2 + ms2)*((-3*g1^2*Lambdax)/5 - 3*g2^2*Lambdax + 
    3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
      Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*Lambdax + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*Lambdax + 
    4*Lambdax^3 - 2*gN^2*Lambdax*QH1p^2 - 2*gN^2*Lambdax*QH2p^2 - 
    2*gN^2*Lambdax*QSp^2 + 3*Lambdax*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    Lambdax*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + 3*Lambdax*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + 
      Yu11^2 + Yu12^2 + Yu20^2 + Yu21^2 + Yu22^2)) + 
  12*gN^2*QQp*QSp*((-2*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 - 
    6*g2^2*MassWB^2 + (g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - 
       mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - 
       mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - 
       ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 
       2*(mu200 + mu211 + mu222)))/5 - 8*gN^2*MassBp^2*QQp^2 + 
    2*gN^2*QQp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TYd02^2 + TYd12^2 + TYd22^2) + 2*(TYu02^2 + TYu12^2 + TYu22^2) + 
    Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
    Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
    Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22) + 
    mq202*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + 
    mq212*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + 
    2*mHd2*(Yd02^2 + Yd12^2 + Yd22^2) + mq222*(Yd02^2 + Yd12^2 + Yd22^2) + 
    2*(Yd02*(md200*Yd02 + md210*Yd12 + md220*Yd22) + 
      Yd12*(md201*Yd02 + md211*Yd12 + md221*Yd22) + 
      Yd22*(md202*Yd02 + md212*Yd12 + md222*Yd22)) + 
    Yu02*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
    Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
    Yu22*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22) + 
    mq202*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 
    mq212*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 
    2*mHu2*(Yu02^2 + Yu12^2 + Yu22^2) + mq222*(Yu02^2 + Yu12^2 + Yu22^2) + 
    2*(Yu02*(mu200*Yu02 + mu210*Yu12 + mu220*Yu22) + 
      Yu12*(mu201*Yu02 + mu211*Yu12 + mu221*Yu22) + 
      Yu22*(mu202*Yu02 + mu212*Yu12 + mu222*Yu22))) + 
  (4*Lambdax^2 + 4*gN^2*QH2p*QSp)*((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 + 
    2*Lambdax^2*mHd2 + 2*Lambdax^2*mHu2 + 2*Lambdax^2*ms2 + 
    (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QH2p^2 + 2*gN^2*QH2p*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    2*TLambdax^2 + 6*(TYu00^2 + TYu01^2 + TYu02^2 + TYu10^2 + TYu11^2 + 
      TYu12^2 + TYu20^2 + TYu21^2 + TYu22^2) + 
    6*mHu2*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + Yu20^2 + 
      Yu21^2 + Yu22^2) + 6*(Yu00*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
      Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
      Yu02*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
      Yu10*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
      Yu11*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
      Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
      Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22) + 
      Yu21*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22) + 
      Yu22*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) + 
    6*(Yu00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
      Yu10*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
      Yu20*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
      Yu01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
      Yu11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
      Yu21*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
      Yu02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22) + 
      Yu12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22) + 
      Yu22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22))) + 
  6*gN^2*QSp*Qup*((-32*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 - 
    (4*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*Qup^2 + 2*gN^2*Qup*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    4*(TYu00^2 + TYu01^2 + TYu02^2) + 4*mHu2*(Yu00^2 + Yu01^2 + Yu02^2) + 
    4*(Yu00*(mq200*Yu00 + mq210*Yu01 + mq220*Yu02) + 
      Yu01*(mq201*Yu00 + mq211*Yu01 + mq221*Yu02) + 
      Yu02*(mq202*Yu00 + mq212*Yu01 + mq222*Yu02)) + 
    2*(Yu00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
      Yu01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
      Yu02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22)) + 
    2*(mu200*(Yu00^2 + Yu01^2 + Yu02^2) + mu210*(Yu00*Yu10 + Yu01*Yu11 + 
        Yu02*Yu12) + mu220*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22))) + 
  6*gN^2*QSp*Qup*((-32*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 - 
    (4*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*Qup^2 + 2*gN^2*Qup*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    4*(TYu10^2 + TYu11^2 + TYu12^2) + 4*mHu2*(Yu10^2 + Yu11^2 + Yu12^2) + 
    4*(Yu10*(mq200*Yu10 + mq210*Yu11 + mq220*Yu12) + 
      Yu11*(mq201*Yu10 + mq211*Yu11 + mq221*Yu12) + 
      Yu12*(mq202*Yu10 + mq212*Yu11 + mq222*Yu12)) + 
    2*(Yu10*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
      Yu11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
      Yu12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22)) + 
    2*(mu201*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + 
      mu211*(Yu10^2 + Yu11^2 + Yu12^2) + mu221*(Yu10*Yu20 + Yu11*Yu21 + 
        Yu12*Yu22))) + 6*gN^2*QSp*Qup*((-32*g1^2*MassB^2)/15 - 
    (32*g3^2*MassG^2)/3 - (4*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - 
       mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - 
       mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - 
       ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 
       2*(mu200 + mu211 + mu222)))/5 - 8*gN^2*MassBp^2*Qup^2 + 
    2*gN^2*Qup*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 
    4*(TYu20^2 + TYu21^2 + TYu22^2) + 4*mHu2*(Yu20^2 + Yu21^2 + Yu22^2) + 
    4*(Yu20*(mq200*Yu20 + mq210*Yu21 + mq220*Yu22) + 
      Yu21*(mq201*Yu20 + mq211*Yu21 + mq221*Yu22) + 
      Yu22*(mq202*Yu20 + mq212*Yu21 + mq222*Yu22)) + 
    2*(Yu20*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
      Yu21*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
      Yu22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) + 
    2*(mu202*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 
      mu212*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 
      mu222*(Yu20^2 + Yu21^2 + Yu22^2))) + 
  8*TLambdax*((6*g1^2*Lambdax*MassB)/5 + 6*g2^2*Lambdax*MassWB + 
    4*gN^2*Lambdax*MassBp*QH1p^2 + 4*gN^2*Lambdax*MassBp*QH2p^2 + 
    4*gN^2*Lambdax*MassBp*QSp^2 + 6*Lambdax*(Kappa00*TKappa00 + 
      Kappa01*TKappa01 + Kappa02*TKappa02 + Kappa10*TKappa10 + 
      Kappa11*TKappa11 + Kappa12*TKappa12 + Kappa20*TKappa20 + 
      Kappa21*TKappa21 + Kappa22*TKappa22) + 
    4*Lambdax*(Lambda1200*TLambda1200 + Lambda1201*TLambda1201 + 
      Lambda1210*TLambda1210 + Lambda1211*TLambda1211) + 
    6*Lambdax*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
      TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
    2*Lambdax*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
      TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 
    6*Lambdax*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + 
      TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) + 
    TLambdax*((-3*g1^2)/5 - 3*g2^2 + 3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + 
        Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + 
        Kappa22^2) + 2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
        Lambda1211^2) + 12*Lambdax^2 - 2*gN^2*QH1p^2 - 2*gN^2*QH2p^2 - 
      2*gN^2*QSp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
        Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + 
      Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2 + 
      3*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + Yu20^2 + 
        Yu21^2 + Yu22^2)))}
