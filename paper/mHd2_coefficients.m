(* Created with the Wolfram Language for Students - Personal Use Only : www.wolfram.com *)
{mHd2, (-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 + 2*Lambdax^2*mHd2 + 
  2*Lambdax^2*mHu2 + 2*Lambdax^2*ms2 - 
  (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
     mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
     mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
     ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
  8*gN^2*MassBp^2*QH1p^2 + 2*gN^2*QH1p*(3*(md200 + md211 + md222)*Qdp + 
    3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
    3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
    2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
    2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
    2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + 
    (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup) + 2*TLambdax^2 + 
  6*(TYd00^2 + TYd01^2 + TYd02^2 + TYd10^2 + TYd11^2 + TYd12^2 + TYd20^2 + 
    TYd21^2 + TYd22^2) + 2*(TYe00^2 + TYe01^2 + TYe02^2 + TYe10^2 + TYe11^2 + 
    TYe12^2 + TYe20^2 + TYe21^2 + TYe22^2) + 
  6*mHd2*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
    Yd21^2 + Yd22^2) + 6*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
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
    Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)), 
 (9*g1^2*g2^2*MassB*MassWB)/5 + (18*g1^2*g2^2*MassWB^2)/5 + 
  87*g2^4*MassWB^2 - 6*Lambdax^2*((Kappa00^2 + Kappa01^2 + Kappa02^2)*
     mDx200 + (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx201 + 
    (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx202 + 
    (Kappa00*Kappa10 + Kappa01*Kappa11 + Kappa02*Kappa12)*mDx210 + 
    (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
    (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx212 + 
    (Kappa00*Kappa20 + Kappa01*Kappa21 + Kappa02*Kappa22)*mDx220 + 
    (Kappa10*Kappa20 + Kappa11*Kappa21 + Kappa12*Kappa22)*mDx221 + 
    (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222) - 
  6*Lambdax^2*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar210 + 
      Kappa02*mDxbar220) + Kappa10*(Kappa10*mDxbar200 + Kappa11*mDxbar210 + 
      Kappa12*mDxbar220) + Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar210 + 
      Kappa22*mDxbar220) + Kappa01*(Kappa00*mDxbar201 + Kappa01*mDxbar211 + 
      Kappa02*mDxbar221) + Kappa11*(Kappa10*mDxbar201 + Kappa11*mDxbar211 + 
      Kappa12*mDxbar221) + Kappa21*(Kappa20*mDxbar201 + Kappa21*mDxbar211 + 
      Kappa22*mDxbar221) + Kappa02*(Kappa00*mDxbar202 + Kappa01*mDxbar212 + 
      Kappa02*mDxbar222) + Kappa12*(Kappa10*mDxbar202 + Kappa11*mDxbar212 + 
      Kappa12*mDxbar222) + Kappa22*(Kappa20*mDxbar202 + Kappa21*mDxbar212 + 
      Kappa22*mDxbar222)) - 4*Lambdax^2*
   (Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I201) + 
    Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I201) + 
    Lambda1201*(Lambda1200*mH1I210 + Lambda1201*mH1I211) + 
    Lambda1211*(Lambda1210*mH1I210 + Lambda1211*mH1I211)) - 
  4*Lambdax^2*((Lambda1200^2 + Lambda1201^2)*mH2I200 + 
    (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I201 + 
    (Lambda1200*Lambda1210 + Lambda1201*Lambda1211)*mH2I210 + 
    (Lambda1210^2 + Lambda1211^2)*mH2I211) - 
  6*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + 
    Kappa20^2 + Kappa21^2 + Kappa22^2)*Lambdax^2*mHd2 - 
  4*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*Lambdax^2*
   mHd2 - 12*Lambdax^4*mHd2 - 6*(Kappa00^2 + Kappa01^2 + Kappa02^2 + 
    Kappa10^2 + Kappa11^2 + Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*
   Lambdax^2*mHu2 - 4*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + 
    Lambda1211^2)*Lambdax^2*mHu2 - 12*Lambdax^4*mHu2 + 
  3*g2^4*(mH1I200 + mH1I211 + mH2I200 + mH2I211 + mHd2 + mHp2 + mHpbar2 + 
    mHu2 + ml200 + ml211 + ml222 + 3*(mq200 + mq211 + mq222)) - 
  12*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + Kappa12^2 + 
    Kappa20^2 + Kappa21^2 + Kappa22^2)*Lambdax^2*ms2 - 
  8*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*Lambdax^2*
   ms2 - 12*Lambdax^4*ms2 + 
  (3*g1^4*(2*(md200 + md211 + md222) + 2*(mDx200 + mDx211 + mDx222) + 
     2*(mDxbar200 + mDxbar211 + mDxbar222) + 6*(me200 + me211 + me222) + 
     3*(mH1I200 + mH1I211) + 3*(mH2I200 + mH2I211) + 3*mHd2 + 3*mHp2 + 
     3*mHpbar2 + 3*mHu2 + 3*(ml200 + ml211 + ml222) + mq200 + mq211 + mq222 + 
     8*(mu200 + mu211 + mu222)))/25 + 12*g2^2*gN^2*MassBp*MassWB*QH1p^2 + 
  24*g2^2*gN^2*MassWB^2*QH1p^2 - 4*gN^2*Lambdax^2*mHd2*QH1p^2 - 
  4*gN^2*Lambdax^2*mHu2*QH1p^2 - 4*gN^2*Lambdax^2*ms2*QH1p^2 + 
  4*gN^2*Lambdax^2*mHd2*QH2p^2 + 4*gN^2*Lambdax^2*mHu2*QH2p^2 + 
  4*gN^2*Lambdax^2*ms2*QH2p^2 + 4*gN^2*Lambdax^2*mHd2*QSp^2 + 
  4*gN^2*Lambdax^2*mHu2*QSp^2 + 4*gN^2*Lambdax^2*ms2*QSp^2 - 
  (24*g1^2*gN^2*QH1p*((md200 + md211 + md222)*Qdp + 
     (mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp - (mDx200 + mDx211 + mDx222)*
      QDxp + (me200 + me211 + me222)*Qep - (mH1I200 + mH1I211)*QH1p - 
     mHd2*QH1p + (mH2I200 + mH2I211)*QH2p + mHu2*QH2p + mHpbar2*QHpbarp - 
     mHp2*QHpp - (ml200 + ml211 + ml222)*QLp + (mq200 + mq211 + mq222)*QQp - 
     2*(mu200 + mu211 + mu222)*Qup))/5 + 
  8*gN^4*QH1p^2*(3*(md200 + md211 + md222)*Qdp^2 + 
    3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp^2 + 
    3*(mDx200 + mDx211 + mDx222)*QDxp^2 + (me200 + me211 + me222)*Qep^2 + 
    2*(mH1I200 + mH1I211)*QH1p^2 + 2*mHd2*QH1p^2 + 
    2*(mH2I200 + mH2I211)*QH2p^2 + 2*mHu2*QH2p^2 + 2*mHpbar2*QHpbarp^2 + 
    2*mHp2*QHpp^2 + 2*(ml200 + ml211 + ml222)*QLp^2 + 
    6*(mq200 + mq211 + mq222)*QQp^2 + ms2*QSp^2 + (msI200 + msI211)*QSp^2 + 
    3*(mu200 + mu211 + mu222)*Qup^2) - 
  6*Lambdax^2*(TKappa00^2 + TKappa01^2 + TKappa02^2 + TKappa10^2 + 
    TKappa11^2 + TKappa12^2 + TKappa20^2 + TKappa21^2 + TKappa22^2) - 
  4*Lambdax^2*(TLambda1200^2 + TLambda1201^2 + TLambda1210^2 + 
    TLambda1211^2) + 4*gN^2*Lambdax*MassBp*QH1p^2*TLambdax - 
  4*gN^2*Lambdax*MassBp*QH2p^2*TLambdax - 4*gN^2*Lambdax*MassBp*QSp^2*
   TLambdax - 12*Lambdax*(Kappa00*TKappa00 + Kappa01*TKappa01 + 
    Kappa02*TKappa02 + Kappa10*TKappa10 + Kappa11*TKappa11 + 
    Kappa12*TKappa12 + Kappa20*TKappa20 + Kappa21*TKappa21 + 
    Kappa22*TKappa22)*TLambdax - 8*Lambdax*(Lambda1200*TLambda1200 + 
    Lambda1201*TLambda1201 + Lambda1210*TLambda1210 + Lambda1211*TLambda1211)*
   TLambdax - 6*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
    Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*TLambdax^2 - 
  4*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*TLambdax^2 - 
  24*Lambdax^2*TLambdax^2 - 4*gN^2*QH1p^2*TLambdax^2 + 
  4*gN^2*QH2p^2*TLambdax^2 + 4*gN^2*QSp^2*TLambdax^2 - 
  (4*g1^2*(TYd00^2 + TYd01^2 + TYd02^2 + TYd10^2 + TYd11^2 + TYd12^2 + 
     TYd20^2 + TYd21^2 + TYd22^2))/5 + 
  32*g3^2*(TYd00^2 + TYd01^2 + TYd02^2 + TYd10^2 + TYd11^2 + TYd12^2 + 
    TYd20^2 + TYd21^2 + TYd22^2) + 12*gN^2*Qdp^2*
   (TYd00^2 + TYd01^2 + TYd02^2 + TYd10^2 + TYd11^2 + TYd12^2 + TYd20^2 + 
    TYd21^2 + TYd22^2) - 12*gN^2*QH1p^2*(TYd00^2 + TYd01^2 + TYd02^2 + 
    TYd10^2 + TYd11^2 + TYd12^2 + TYd20^2 + TYd21^2 + TYd22^2) + 
  12*gN^2*QQp^2*(TYd00^2 + TYd01^2 + TYd02^2 + TYd10^2 + TYd11^2 + TYd12^2 + 
    TYd20^2 + TYd21^2 + TYd22^2) + 
  (12*g1^2*(TYe00^2 + TYe01^2 + TYe02^2 + TYe10^2 + TYe11^2 + TYe12^2 + 
     TYe20^2 + TYe21^2 + TYe22^2))/5 + 
  4*gN^2*Qep^2*(TYe00^2 + TYe01^2 + TYe02^2 + TYe10^2 + TYe11^2 + TYe12^2 + 
    TYe20^2 + TYe21^2 + TYe22^2) - 4*gN^2*QH1p^2*
   (TYe00^2 + TYe01^2 + TYe02^2 + TYe10^2 + TYe11^2 + TYe12^2 + TYe20^2 + 
    TYe21^2 + TYe22^2) + 4*gN^2*QLp^2*(TYe00^2 + TYe01^2 + TYe02^2 + 
    TYe10^2 + TYe11^2 + TYe12^2 + TYe20^2 + TYe21^2 + TYe22^2) - 
  6*Lambdax^2*(TYu00^2 + TYu01^2 + TYu02^2 + TYu10^2 + TYu11^2 + TYu12^2 + 
    TYu20^2 + TYu21^2 + TYu22^2) + 
  (4*g1^2*MassB*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
     TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22))/5 - 
  64*g3^2*MassG*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
    TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) - 
  12*gN^2*MassBp*Qdp^2*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
    TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
  12*gN^2*MassBp*QH1p^2*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
    TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) - 
  12*gN^2*MassBp*QQp^2*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
    TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
  64*g3^2*MassG^2*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
    Yd20^2 + Yd21^2 + Yd22^2) - 
  (4*g1^2*mHd2*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
     Yd20^2 + Yd21^2 + Yd22^2))/5 + 32*g3^2*mHd2*(Yd00^2 + Yd01^2 + Yd02^2 + 
    Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
  12*gN^2*mHd2*Qdp^2*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
    Yd20^2 + Yd21^2 + Yd22^2) - 12*gN^2*mHd2*QH1p^2*
   (Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + 
    Yd22^2) + 12*gN^2*mHd2*QQp^2*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
    Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) - 
  (4*g1^2*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
     Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
     Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
     Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
     Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
     Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
     Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22) + 
     Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22) + 
     Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22)))/5 + 
  32*g3^2*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
    Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
    Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
    Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
    Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
    Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
    Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22) + 
    Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22) + 
    Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
  12*gN^2*Qdp^2*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
    Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
    Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
    Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
    Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
    Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
    Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22) + 
    Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22) + 
    Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22)) - 
  12*gN^2*QH1p^2*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
    Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
    Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
    Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
    Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
    Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
    Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22) + 
    Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22) + 
    Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
  12*gN^2*QQp^2*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
    Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
    Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
    Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
    Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
    Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
    Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22) + 
    Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22) + 
    Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22)) - 
  (4*g1^2*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
     Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
     Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
     Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
     Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
     Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
     Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + 
     Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + 
     Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)))/5 + 
  32*g3^2*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
    Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
    Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
    Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
    Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
    Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
    Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + 
    Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + 
    Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) + 
  12*gN^2*Qdp^2*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
    Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
    Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
    Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
    Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
    Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
    Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + 
    Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + 
    Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) - 
  12*gN^2*QH1p^2*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
    Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
    Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
    Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
    Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
    Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
    Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + 
    Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + 
    Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) + 
  12*gN^2*QQp^2*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
    Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
    Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
    Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
    Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
    Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
    Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + 
    Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + 
    Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) - 
  36*(Yd00*(TYd00*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02) + 
      TYd10*(TYd10*Yd00 + TYd11*Yd01 + TYd12*Yd02) + 
      TYd20*(TYd20*Yd00 + TYd21*Yd01 + TYd22*Yd02)) + 
    Yd01*(TYd01*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02) + 
      TYd11*(TYd10*Yd00 + TYd11*Yd01 + TYd12*Yd02) + 
      TYd21*(TYd20*Yd00 + TYd21*Yd01 + TYd22*Yd02)) + 
    Yd02*(TYd02*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02) + 
      TYd12*(TYd10*Yd00 + TYd11*Yd01 + TYd12*Yd02) + 
      TYd22*(TYd20*Yd00 + TYd21*Yd01 + TYd22*Yd02)) + 
    Yd10*(TYd00*(TYd00*Yd10 + TYd01*Yd11 + TYd02*Yd12) + 
      TYd10*(TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12) + 
      TYd20*(TYd20*Yd10 + TYd21*Yd11 + TYd22*Yd12)) + 
    Yd11*(TYd01*(TYd00*Yd10 + TYd01*Yd11 + TYd02*Yd12) + 
      TYd11*(TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12) + 
      TYd21*(TYd20*Yd10 + TYd21*Yd11 + TYd22*Yd12)) + 
    Yd12*(TYd02*(TYd00*Yd10 + TYd01*Yd11 + TYd02*Yd12) + 
      TYd12*(TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12) + 
      TYd22*(TYd20*Yd10 + TYd21*Yd11 + TYd22*Yd12)) + 
    Yd20*(TYd00*(TYd00*Yd20 + TYd01*Yd21 + TYd02*Yd22) + 
      TYd10*(TYd10*Yd20 + TYd11*Yd21 + TYd12*Yd22) + 
      TYd20*(TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)) + 
    Yd21*(TYd01*(TYd00*Yd20 + TYd01*Yd21 + TYd02*Yd22) + 
      TYd11*(TYd10*Yd20 + TYd11*Yd21 + TYd12*Yd22) + 
      TYd21*(TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)) + 
    Yd22*(TYd02*(TYd00*Yd20 + TYd01*Yd21 + TYd02*Yd22) + 
      TYd12*(TYd10*Yd20 + TYd11*Yd21 + TYd12*Yd22) + 
      TYd22*(TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22))) - 
  6*(Yd00*(TYu00*(TYu00*Yd00 + TYu01*Yd01 + TYu02*Yd02) + 
      TYu10*(TYu10*Yd00 + TYu11*Yd01 + TYu12*Yd02) + 
      TYu20*(TYu20*Yd00 + TYu21*Yd01 + TYu22*Yd02)) + 
    Yd01*(TYu01*(TYu00*Yd00 + TYu01*Yd01 + TYu02*Yd02) + 
      TYu11*(TYu10*Yd00 + TYu11*Yd01 + TYu12*Yd02) + 
      TYu21*(TYu20*Yd00 + TYu21*Yd01 + TYu22*Yd02)) + 
    Yd02*(TYu02*(TYu00*Yd00 + TYu01*Yd01 + TYu02*Yd02) + 
      TYu12*(TYu10*Yd00 + TYu11*Yd01 + TYu12*Yd02) + 
      TYu22*(TYu20*Yd00 + TYu21*Yd01 + TYu22*Yd02)) + 
    Yd10*(TYu00*(TYu00*Yd10 + TYu01*Yd11 + TYu02*Yd12) + 
      TYu10*(TYu10*Yd10 + TYu11*Yd11 + TYu12*Yd12) + 
      TYu20*(TYu20*Yd10 + TYu21*Yd11 + TYu22*Yd12)) + 
    Yd11*(TYu01*(TYu00*Yd10 + TYu01*Yd11 + TYu02*Yd12) + 
      TYu11*(TYu10*Yd10 + TYu11*Yd11 + TYu12*Yd12) + 
      TYu21*(TYu20*Yd10 + TYu21*Yd11 + TYu22*Yd12)) + 
    Yd12*(TYu02*(TYu00*Yd10 + TYu01*Yd11 + TYu02*Yd12) + 
      TYu12*(TYu10*Yd10 + TYu11*Yd11 + TYu12*Yd12) + 
      TYu22*(TYu20*Yd10 + TYu21*Yd11 + TYu22*Yd12)) + 
    Yd20*(TYu00*(TYu00*Yd20 + TYu01*Yd21 + TYu02*Yd22) + 
      TYu10*(TYu10*Yd20 + TYu11*Yd21 + TYu12*Yd22) + 
      TYu20*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22)) + 
    Yd21*(TYu01*(TYu00*Yd20 + TYu01*Yd21 + TYu02*Yd22) + 
      TYu11*(TYu10*Yd20 + TYu11*Yd21 + TYu12*Yd22) + 
      TYu21*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22)) + 
    Yd22*(TYu02*(TYu00*Yd20 + TYu01*Yd21 + TYu02*Yd22) + 
      TYu12*(TYu10*Yd20 + TYu11*Yd21 + TYu12*Yd22) + 
      TYu22*(TYu20*Yd20 + TYu21*Yd21 + TYu22*Yd22))) - 
  36*(TYd00*(TYd00*(Yd00^2 + Yd01^2 + Yd02^2) + 
      TYd10*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      TYd20*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    TYd01*(TYd01*(Yd00^2 + Yd01^2 + Yd02^2) + 
      TYd11*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      TYd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    TYd02*(TYd02*(Yd00^2 + Yd01^2 + Yd02^2) + 
      TYd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      TYd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    TYd10*(TYd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      TYd10*(Yd10^2 + Yd11^2 + Yd12^2) + TYd20*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + TYd11*(TYd01*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      TYd11*(Yd10^2 + Yd11^2 + Yd12^2) + TYd21*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + TYd12*(TYd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      TYd12*(Yd10^2 + Yd11^2 + Yd12^2) + TYd22*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + TYd20*(TYd00*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      TYd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      TYd20*(Yd20^2 + Yd21^2 + Yd22^2)) + 
    TYd21*(TYd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      TYd11*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      TYd21*(Yd20^2 + Yd21^2 + Yd22^2)) + 
    TYd22*(TYd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      TYd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      TYd22*(Yd20^2 + Yd21^2 + Yd22^2))) - 
  36*mHd2*(Yd00*(Yd00*(Yd00^2 + Yd01^2 + Yd02^2) + 
      Yd10*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd20*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    Yd01*(Yd01*(Yd00^2 + Yd01^2 + Yd02^2) + 
      Yd11*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    Yd02*(Yd02*(Yd00^2 + Yd01^2 + Yd02^2) + 
      Yd12*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    Yd10*(Yd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd10*(Yd10^2 + Yd11^2 + Yd12^2) + Yd20*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + Yd11*(Yd01*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd11*(Yd10^2 + Yd11^2 + Yd12^2) + Yd21*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + Yd12*(Yd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd12*(Yd10^2 + Yd11^2 + Yd12^2) + Yd22*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + Yd20*(Yd00*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      Yd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      Yd20*(Yd20^2 + Yd21^2 + Yd22^2)) + 
    Yd21*(Yd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      Yd11*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      Yd21*(Yd20^2 + Yd21^2 + Yd22^2)) + 
    Yd22*(Yd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      Yd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      Yd22*(Yd20^2 + Yd21^2 + Yd22^2))) - 
  36*(Yd00*(Yd00*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
        Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
        Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22)) + 
      Yd10*(Yd10*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
        Yd11*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
        Yd12*(md200*Yd02 + md201*Yd12 + md202*Yd22)) + 
      Yd20*(Yd20*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
        Yd21*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
        Yd22*(md200*Yd02 + md201*Yd12 + md202*Yd22))) + 
    Yd01*(Yd01*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
        Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
        Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22)) + 
      Yd11*(Yd10*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
        Yd11*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
        Yd12*(md200*Yd02 + md201*Yd12 + md202*Yd22)) + 
      Yd21*(Yd20*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
        Yd21*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
        Yd22*(md200*Yd02 + md201*Yd12 + md202*Yd22))) + 
    Yd02*(Yd02*(Yd00*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
        Yd01*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
        Yd02*(md200*Yd02 + md201*Yd12 + md202*Yd22)) + 
      Yd12*(Yd10*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
        Yd11*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
        Yd12*(md200*Yd02 + md201*Yd12 + md202*Yd22)) + 
      Yd22*(Yd20*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
        Yd21*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
        Yd22*(md200*Yd02 + md201*Yd12 + md202*Yd22))) + 
    Yd10*(Yd00*(Yd00*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
        Yd01*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
        Yd02*(md210*Yd02 + md211*Yd12 + md212*Yd22)) + 
      Yd10*(Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
        Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
        Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22)) + 
      Yd20*(Yd20*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
        Yd21*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
        Yd22*(md210*Yd02 + md211*Yd12 + md212*Yd22))) + 
    Yd11*(Yd01*(Yd00*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
        Yd01*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
        Yd02*(md210*Yd02 + md211*Yd12 + md212*Yd22)) + 
      Yd11*(Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
        Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
        Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22)) + 
      Yd21*(Yd20*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
        Yd21*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
        Yd22*(md210*Yd02 + md211*Yd12 + md212*Yd22))) + 
    Yd12*(Yd02*(Yd00*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
        Yd01*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
        Yd02*(md210*Yd02 + md211*Yd12 + md212*Yd22)) + 
      Yd12*(Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
        Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
        Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22)) + 
      Yd22*(Yd20*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
        Yd21*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
        Yd22*(md210*Yd02 + md211*Yd12 + md212*Yd22))) + 
    Yd20*(Yd00*(Yd00*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
        Yd01*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
        Yd02*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
      Yd10*(Yd10*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
        Yd11*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
        Yd12*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
      Yd20*(Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
        Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
        Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22))) + 
    Yd21*(Yd01*(Yd00*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
        Yd01*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
        Yd02*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
      Yd11*(Yd10*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
        Yd11*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
        Yd12*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
      Yd21*(Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
        Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
        Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22))) + 
    Yd22*(Yd02*(Yd00*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
        Yd01*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
        Yd02*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
      Yd12*(Yd10*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
        Yd11*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
        Yd12*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
      Yd22*(Yd20*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
        Yd21*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
        Yd22*(md220*Yd02 + md221*Yd12 + md222*Yd22)))) - 
  36*(Yd00*(Yd00*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22)) + 
      Yd01*(Yd01*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd11*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd21*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22)) + 
      Yd02*(Yd02*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd12*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd22*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))) + 
    Yd10*(Yd10*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22)) + 
      Yd11*(Yd01*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd11*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd21*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22)) + 
      Yd12*(Yd02*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd12*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd22*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))) + 
    Yd20*(Yd20*(Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22)) + 
      Yd21*(Yd01*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd11*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd21*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22)) + 
      Yd22*(Yd02*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd12*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd22*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))) + 
    Yd01*(Yd00*(Yd00*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd10*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd20*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22)) + 
      Yd01*(Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22)) + 
      Yd02*(Yd02*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd12*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd22*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))) + 
    Yd11*(Yd10*(Yd00*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd10*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd20*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22)) + 
      Yd11*(Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22)) + 
      Yd12*(Yd02*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd12*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd22*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))) + 
    Yd21*(Yd20*(Yd00*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd10*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd20*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22)) + 
      Yd21*(Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22)) + 
      Yd22*(Yd02*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd12*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd22*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))) + 
    Yd02*(Yd00*(Yd00*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd10*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd20*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) + 
      Yd01*(Yd01*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd11*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd21*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) + 
      Yd02*(Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))) + 
    Yd12*(Yd10*(Yd00*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd10*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd20*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) + 
      Yd11*(Yd01*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd11*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd21*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) + 
      Yd12*(Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))) + 
    Yd22*(Yd20*(Yd00*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd10*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd20*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) + 
      Yd21*(Yd01*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd11*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd21*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)) + 
      Yd22*(Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22)))) - 
  (12*g1^2*MassB*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
     TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22))/5 - 
  4*gN^2*MassBp*Qep^2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
    TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 
  4*gN^2*MassBp*QH1p^2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
    TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) - 
  4*gN^2*MassBp*QLp^2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
    TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 
  (12*g1^2*mHd2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + 
     Ye20^2 + Ye21^2 + Ye22^2))/5 + 4*gN^2*mHd2*Qep^2*
   (Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + 
    Ye22^2) - 4*gN^2*mHd2*QH1p^2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + 
    Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
  4*gN^2*mHd2*QLp^2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + 
    Ye20^2 + Ye21^2 + Ye22^2) + 
  (12*g1^2*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
     Ye10*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
     Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
     Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
     Ye11*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
     Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
     Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22) + 
     Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22) + 
     Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)))/5 + 
  4*gN^2*Qep^2*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
    Ye10*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
    Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
    Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
    Ye11*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
    Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
    Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22) + 
    Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22) + 
    Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)) - 
  4*gN^2*QH1p^2*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
    Ye10*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
    Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
    Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
    Ye11*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
    Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
    Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22) + 
    Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22) + 
    Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
  4*gN^2*QLp^2*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
    Ye10*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
    Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
    Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
    Ye11*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
    Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
    Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22) + 
    Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22) + 
    Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
  (12*g1^2*(Ye00*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
     Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
     Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
     Ye10*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
     Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
     Ye12*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
     Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22) + 
     Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + 
     Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)))/5 + 
  4*gN^2*Qep^2*(Ye00*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
    Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
    Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
    Ye10*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
    Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
    Ye12*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
    Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22) + 
    Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + 
    Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)) - 
  4*gN^2*QH1p^2*(Ye00*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
    Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
    Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
    Ye10*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
    Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
    Ye12*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
    Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22) + 
    Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + 
    Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)) + 
  4*gN^2*QLp^2*(Ye00*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
    Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
    Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
    Ye10*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
    Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
    Ye12*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
    Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22) + 
    Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + 
    Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)) + 
  (g1^2*MassB*(891*g1^2*MassB + 90*g2^2*MassB + 45*g2^2*MassWB - 
     360*gN^2*MassB*Qdp*QH1p - 180*gN^2*MassBp*Qdp*QH1p - 
     360*gN^2*MassB*QDxbarp*QH1p - 180*gN^2*MassBp*QDxbarp*QH1p + 
     360*gN^2*MassB*QDxp*QH1p + 180*gN^2*MassBp*QDxp*QH1p - 
     360*gN^2*MassB*Qep*QH1p - 180*gN^2*MassBp*Qep*QH1p + 
     480*gN^2*MassB*QH1p^2 + 240*gN^2*MassBp*QH1p^2 - 
     360*gN^2*MassB*QH1p*QH2p - 180*gN^2*MassBp*QH1p*QH2p - 
     120*gN^2*MassB*QH1p*QHpbarp - 60*gN^2*MassBp*QH1p*QHpbarp + 
     120*gN^2*MassB*QH1p*QHpp + 60*gN^2*MassBp*QH1p*QHpp + 
     360*gN^2*MassB*QH1p*QLp + 180*gN^2*MassBp*QH1p*QLp - 
     360*gN^2*MassB*QH1p*QQp - 180*gN^2*MassBp*QH1p*QQp + 
     720*gN^2*MassB*QH1p*Qup + 360*gN^2*MassBp*QH1p*Qup + 
     20*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
       TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) - 
     40*MassB*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
       Yd21^2 + Yd22^2) - 60*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + 
       TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + 
       TYe22*Ye22) + 120*MassB*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
       Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2)))/25 + 
  (4*gN^2*MassBp*(-9*g1^2*MassB*Qdp*QH1p - 18*g1^2*MassBp*Qdp*QH1p - 
     9*g1^2*MassB*QDxbarp*QH1p - 18*g1^2*MassBp*QDxbarp*QH1p + 
     9*g1^2*MassB*QDxp*QH1p + 18*g1^2*MassBp*QDxp*QH1p - 
     9*g1^2*MassB*Qep*QH1p - 18*g1^2*MassBp*Qep*QH1p + 12*g1^2*MassB*QH1p^2 + 
     24*g1^2*MassBp*QH1p^2 + 30*g2^2*MassBp*QH1p^2 + 15*g2^2*MassWB*QH1p^2 + 
     270*gN^2*MassBp*Qdp^2*QH1p^2 + 270*gN^2*MassBp*QDxbarp^2*QH1p^2 + 
     270*gN^2*MassBp*QDxp^2*QH1p^2 + 90*gN^2*MassBp*Qep^2*QH1p^2 + 
     240*gN^2*MassBp*QH1p^4 - 9*g1^2*MassB*QH1p*QH2p - 
     18*g1^2*MassBp*QH1p*QH2p + 180*gN^2*MassBp*QH1p^2*QH2p^2 - 
     3*g1^2*MassB*QH1p*QHpbarp - 6*g1^2*MassBp*QH1p*QHpbarp + 
     60*gN^2*MassBp*QH1p^2*QHpbarp^2 + 3*g1^2*MassB*QH1p*QHpp + 
     6*g1^2*MassBp*QH1p*QHpp + 60*gN^2*MassBp*QH1p^2*QHpp^2 + 
     9*g1^2*MassB*QH1p*QLp + 18*g1^2*MassBp*QH1p*QLp + 
     180*gN^2*MassBp*QH1p^2*QLp^2 - 9*g1^2*MassB*QH1p*QQp - 
     18*g1^2*MassBp*QH1p*QQp + 540*gN^2*MassBp*QH1p^2*QQp^2 + 
     90*gN^2*MassBp*QH1p^2*QSp^2 + 18*g1^2*MassB*QH1p*Qup + 
     36*g1^2*MassBp*QH1p*Qup + 270*gN^2*MassBp*QH1p^2*Qup^2 - 
     5*Lambdax*(QH1p^2 - QH2p^2 - QSp^2)*(2*Lambdax*MassBp - TLambdax) - 
     15*Qdp^2*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
       TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
     15*QH1p^2*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
       TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) - 
     15*QQp^2*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + 
       TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
     30*MassBp*(Qdp^2 - QH1p^2 + QQp^2)*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
       Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) - 
     5*Qep^2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
       TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 
     5*QH1p^2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
       TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) - 
     5*QLp^2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
       TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22) + 
     10*MassBp*Qep^2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + 
       Ye20^2 + Ye21^2 + Ye22^2) - 10*MassBp*QH1p^2*(Ye00^2 + Ye01^2 + 
       Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
     10*MassBp*QLp^2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + 
       Ye20^2 + Ye21^2 + Ye22^2)))/5 - 
  12*(Ye00*(TYe00*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02) + 
      TYe10*(TYe10*Ye00 + TYe11*Ye01 + TYe12*Ye02) + 
      TYe20*(TYe20*Ye00 + TYe21*Ye01 + TYe22*Ye02)) + 
    Ye01*(TYe01*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02) + 
      TYe11*(TYe10*Ye00 + TYe11*Ye01 + TYe12*Ye02) + 
      TYe21*(TYe20*Ye00 + TYe21*Ye01 + TYe22*Ye02)) + 
    Ye02*(TYe02*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02) + 
      TYe12*(TYe10*Ye00 + TYe11*Ye01 + TYe12*Ye02) + 
      TYe22*(TYe20*Ye00 + TYe21*Ye01 + TYe22*Ye02)) + 
    Ye10*(TYe00*(TYe00*Ye10 + TYe01*Ye11 + TYe02*Ye12) + 
      TYe10*(TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12) + 
      TYe20*(TYe20*Ye10 + TYe21*Ye11 + TYe22*Ye12)) + 
    Ye11*(TYe01*(TYe00*Ye10 + TYe01*Ye11 + TYe02*Ye12) + 
      TYe11*(TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12) + 
      TYe21*(TYe20*Ye10 + TYe21*Ye11 + TYe22*Ye12)) + 
    Ye12*(TYe02*(TYe00*Ye10 + TYe01*Ye11 + TYe02*Ye12) + 
      TYe12*(TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12) + 
      TYe22*(TYe20*Ye10 + TYe21*Ye11 + TYe22*Ye12)) + 
    Ye20*(TYe00*(TYe00*Ye20 + TYe01*Ye21 + TYe02*Ye22) + 
      TYe10*(TYe10*Ye20 + TYe11*Ye21 + TYe12*Ye22) + 
      TYe20*(TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Ye21*(TYe01*(TYe00*Ye20 + TYe01*Ye21 + TYe02*Ye22) + 
      TYe11*(TYe10*Ye20 + TYe11*Ye21 + TYe12*Ye22) + 
      TYe21*(TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Ye22*(TYe02*(TYe00*Ye20 + TYe01*Ye21 + TYe02*Ye22) + 
      TYe12*(TYe10*Ye20 + TYe11*Ye21 + TYe12*Ye22) + 
      TYe22*(TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22))) - 
  12*(TYe00*(TYe00*(Ye00^2 + Ye01^2 + Ye02^2) + 
      TYe10*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      TYe20*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22)) + 
    TYe01*(TYe01*(Ye00^2 + Ye01^2 + Ye02^2) + 
      TYe11*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      TYe21*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22)) + 
    TYe02*(TYe02*(Ye00^2 + Ye01^2 + Ye02^2) + 
      TYe12*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      TYe22*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22)) + 
    TYe10*(TYe00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      TYe10*(Ye10^2 + Ye11^2 + Ye12^2) + TYe20*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22)) + TYe11*(TYe01*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      TYe11*(Ye10^2 + Ye11^2 + Ye12^2) + TYe21*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22)) + TYe12*(TYe02*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      TYe12*(Ye10^2 + Ye11^2 + Ye12^2) + TYe22*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22)) + TYe20*(TYe00*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      TYe10*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      TYe20*(Ye20^2 + Ye21^2 + Ye22^2)) + 
    TYe21*(TYe01*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      TYe11*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      TYe21*(Ye20^2 + Ye21^2 + Ye22^2)) + 
    TYe22*(TYe02*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      TYe12*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      TYe22*(Ye20^2 + Ye21^2 + Ye22^2))) - 
  12*mHd2*(Ye00*(Ye00*(Ye00^2 + Ye01^2 + Ye02^2) + 
      Ye10*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      Ye20*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22)) + 
    Ye01*(Ye01*(Ye00^2 + Ye01^2 + Ye02^2) + 
      Ye11*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      Ye21*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22)) + 
    Ye02*(Ye02*(Ye00^2 + Ye01^2 + Ye02^2) + 
      Ye12*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      Ye22*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22)) + 
    Ye10*(Ye00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      Ye10*(Ye10^2 + Ye11^2 + Ye12^2) + Ye20*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22)) + Ye11*(Ye01*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      Ye11*(Ye10^2 + Ye11^2 + Ye12^2) + Ye21*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22)) + Ye12*(Ye02*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      Ye12*(Ye10^2 + Ye11^2 + Ye12^2) + Ye22*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22)) + Ye20*(Ye00*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      Ye10*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      Ye20*(Ye20^2 + Ye21^2 + Ye22^2)) + 
    Ye21*(Ye01*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      Ye11*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      Ye21*(Ye20^2 + Ye21^2 + Ye22^2)) + 
    Ye22*(Ye02*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      Ye12*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      Ye22*(Ye20^2 + Ye21^2 + Ye22^2))) - 
  12*(Ye00*(Ye00*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
        Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
        Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22)) + 
      Ye10*(Ye10*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
        Ye11*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
        Ye12*(me200*Ye02 + me201*Ye12 + me202*Ye22)) + 
      Ye20*(Ye20*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
        Ye21*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
        Ye22*(me200*Ye02 + me201*Ye12 + me202*Ye22))) + 
    Ye01*(Ye01*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
        Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
        Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22)) + 
      Ye11*(Ye10*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
        Ye11*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
        Ye12*(me200*Ye02 + me201*Ye12 + me202*Ye22)) + 
      Ye21*(Ye20*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
        Ye21*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
        Ye22*(me200*Ye02 + me201*Ye12 + me202*Ye22))) + 
    Ye02*(Ye02*(Ye00*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
        Ye01*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
        Ye02*(me200*Ye02 + me201*Ye12 + me202*Ye22)) + 
      Ye12*(Ye10*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
        Ye11*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
        Ye12*(me200*Ye02 + me201*Ye12 + me202*Ye22)) + 
      Ye22*(Ye20*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
        Ye21*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
        Ye22*(me200*Ye02 + me201*Ye12 + me202*Ye22))) + 
    Ye10*(Ye00*(Ye00*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
        Ye01*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
        Ye02*(me210*Ye02 + me211*Ye12 + me212*Ye22)) + 
      Ye10*(Ye10*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
        Ye11*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
        Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22)) + 
      Ye20*(Ye20*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
        Ye21*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
        Ye22*(me210*Ye02 + me211*Ye12 + me212*Ye22))) + 
    Ye11*(Ye01*(Ye00*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
        Ye01*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
        Ye02*(me210*Ye02 + me211*Ye12 + me212*Ye22)) + 
      Ye11*(Ye10*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
        Ye11*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
        Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22)) + 
      Ye21*(Ye20*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
        Ye21*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
        Ye22*(me210*Ye02 + me211*Ye12 + me212*Ye22))) + 
    Ye12*(Ye02*(Ye00*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
        Ye01*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
        Ye02*(me210*Ye02 + me211*Ye12 + me212*Ye22)) + 
      Ye12*(Ye10*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
        Ye11*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
        Ye12*(me210*Ye02 + me211*Ye12 + me212*Ye22)) + 
      Ye22*(Ye20*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
        Ye21*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
        Ye22*(me210*Ye02 + me211*Ye12 + me212*Ye22))) + 
    Ye20*(Ye00*(Ye00*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
        Ye01*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
        Ye02*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
      Ye10*(Ye10*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
        Ye11*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
        Ye12*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
      Ye20*(Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
        Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
        Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22))) + 
    Ye21*(Ye01*(Ye00*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
        Ye01*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
        Ye02*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
      Ye11*(Ye10*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
        Ye11*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
        Ye12*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
      Ye21*(Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
        Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
        Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22))) + 
    Ye22*(Ye02*(Ye00*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
        Ye01*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
        Ye02*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
      Ye12*(Ye10*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
        Ye11*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
        Ye12*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
      Ye22*(Ye20*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
        Ye21*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
        Ye22*(me220*Ye02 + me221*Ye12 + me222*Ye22)))) - 
  12*(Ye00*(Ye00*(Ye00*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
        Ye10*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
        Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22)) + 
      Ye01*(Ye01*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
        Ye11*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
        Ye21*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22)) + 
      Ye02*(Ye02*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
        Ye12*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
        Ye22*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22))) + 
    Ye10*(Ye10*(Ye00*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
        Ye10*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
        Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22)) + 
      Ye11*(Ye01*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
        Ye11*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
        Ye21*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22)) + 
      Ye12*(Ye02*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
        Ye12*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
        Ye22*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22))) + 
    Ye20*(Ye20*(Ye00*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
        Ye10*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
        Ye20*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22)) + 
      Ye21*(Ye01*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
        Ye11*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
        Ye21*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22)) + 
      Ye22*(Ye02*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
        Ye12*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
        Ye22*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22))) + 
    Ye01*(Ye00*(Ye00*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
        Ye10*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
        Ye20*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22)) + 
      Ye01*(Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
        Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
        Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22)) + 
      Ye02*(Ye02*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
        Ye12*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
        Ye22*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22))) + 
    Ye11*(Ye10*(Ye00*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
        Ye10*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
        Ye20*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22)) + 
      Ye11*(Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
        Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
        Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22)) + 
      Ye12*(Ye02*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
        Ye12*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
        Ye22*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22))) + 
    Ye21*(Ye20*(Ye00*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
        Ye10*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
        Ye20*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22)) + 
      Ye21*(Ye01*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
        Ye11*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
        Ye21*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22)) + 
      Ye22*(Ye02*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
        Ye12*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
        Ye22*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22))) + 
    Ye02*(Ye00*(Ye00*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
        Ye10*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
        Ye20*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)) + 
      Ye01*(Ye01*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
        Ye11*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
        Ye21*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)) + 
      Ye02*(Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
        Ye12*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
        Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22))) + 
    Ye12*(Ye10*(Ye00*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
        Ye10*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
        Ye20*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)) + 
      Ye11*(Ye01*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
        Ye11*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
        Ye21*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)) + 
      Ye12*(Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
        Ye12*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
        Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22))) + 
    Ye22*(Ye20*(Ye00*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
        Ye10*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
        Ye20*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)) + 
      Ye21*(Ye01*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
        Ye11*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
        Ye21*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)) + 
      Ye22*(Ye02*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
        Ye12*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
        Ye22*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22)))) - 
  12*Lambdax*TLambdax*(TYu00*Yu00 + TYu01*Yu01 + TYu02*Yu02 + TYu10*Yu10 + 
    TYu11*Yu11 + TYu12*Yu12 + TYu20*Yu20 + TYu21*Yu21 + TYu22*Yu22) - 
  6*Lambdax^2*mHd2*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + 
    Yu20^2 + Yu21^2 + Yu22^2) - 12*Lambdax^2*mHu2*
   (Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + Yu20^2 + Yu21^2 + 
    Yu22^2) - 6*Lambdax^2*ms2*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + 
    Yu12^2 + Yu20^2 + Yu21^2 + Yu22^2) - 
  6*TLambdax^2*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + 
    Yu20^2 + Yu21^2 + Yu22^2) - 6*Lambdax^2*
   (Yu00*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
    Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
    Yu02*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
    Yu10*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
    Yu11*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
    Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
    Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22) + 
    Yu21*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22) + 
    Yu22*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) - 
  6*Lambdax^2*(Yu00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
    Yu10*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
    Yu20*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
    Yu01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
    Yu11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
    Yu21*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
    Yu02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22) + 
    Yu12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22) + 
    Yu22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) - 
  6*(Yu00*((Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))*Yu00 + 
      (Yd01*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd11*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd21*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))*Yu01 + 
      (Yd02*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd12*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd22*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))*Yu02) + 
    Yu01*((Yd00*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd10*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd20*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))*Yu00 + 
      (Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))*Yu01 + 
      (Yd02*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd12*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd22*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))*Yu02) + 
    Yu02*((Yd00*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd10*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd20*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))*Yu00 + 
      (Yd01*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd11*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd21*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))*Yu01 + 
      (Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))*Yu02) + 
    Yu10*((Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))*Yu10 + 
      (Yd01*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd11*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd21*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))*Yu11 + 
      (Yd02*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd12*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd22*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))*Yu12) + 
    Yu11*((Yd00*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd10*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd20*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))*Yu10 + 
      (Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))*Yu11 + 
      (Yd02*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd12*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd22*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))*Yu12) + 
    Yu12*((Yd00*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd10*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd20*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))*Yu10 + 
      (Yd01*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd11*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd21*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))*Yu11 + 
      (Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))*Yu12) + 
    Yu20*((Yd00*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd10*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd20*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))*Yu20 + 
      (Yd01*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd11*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd21*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))*Yu21 + 
      (Yd02*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
        Yd12*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
        Yd22*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22))*Yu22) + 
    Yu21*((Yd00*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd10*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd20*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))*Yu20 + 
      (Yd01*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd11*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd21*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))*Yu21 + 
      (Yd02*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
        Yd12*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
        Yd22*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22))*Yu22) + 
    Yu22*((Yd00*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd10*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd20*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))*Yu20 + 
      (Yd01*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd11*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd21*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))*Yu21 + 
      (Yd02*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
        Yd12*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
        Yd22*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22))*Yu22)) - 
  6*(Yu00*(TYd00*(TYd00*Yu00 + TYd01*Yu01 + TYd02*Yu02) + 
      TYd10*(TYd10*Yu00 + TYd11*Yu01 + TYd12*Yu02) + 
      TYd20*(TYd20*Yu00 + TYd21*Yu01 + TYd22*Yu02)) + 
    Yu01*(TYd01*(TYd00*Yu00 + TYd01*Yu01 + TYd02*Yu02) + 
      TYd11*(TYd10*Yu00 + TYd11*Yu01 + TYd12*Yu02) + 
      TYd21*(TYd20*Yu00 + TYd21*Yu01 + TYd22*Yu02)) + 
    Yu02*(TYd02*(TYd00*Yu00 + TYd01*Yu01 + TYd02*Yu02) + 
      TYd12*(TYd10*Yu00 + TYd11*Yu01 + TYd12*Yu02) + 
      TYd22*(TYd20*Yu00 + TYd21*Yu01 + TYd22*Yu02)) + 
    Yu10*(TYd00*(TYd00*Yu10 + TYd01*Yu11 + TYd02*Yu12) + 
      TYd10*(TYd10*Yu10 + TYd11*Yu11 + TYd12*Yu12) + 
      TYd20*(TYd20*Yu10 + TYd21*Yu11 + TYd22*Yu12)) + 
    Yu11*(TYd01*(TYd00*Yu10 + TYd01*Yu11 + TYd02*Yu12) + 
      TYd11*(TYd10*Yu10 + TYd11*Yu11 + TYd12*Yu12) + 
      TYd21*(TYd20*Yu10 + TYd21*Yu11 + TYd22*Yu12)) + 
    Yu12*(TYd02*(TYd00*Yu10 + TYd01*Yu11 + TYd02*Yu12) + 
      TYd12*(TYd10*Yu10 + TYd11*Yu11 + TYd12*Yu12) + 
      TYd22*(TYd20*Yu10 + TYd21*Yu11 + TYd22*Yu12)) + 
    Yu20*(TYd00*(TYd00*Yu20 + TYd01*Yu21 + TYd02*Yu22) + 
      TYd10*(TYd10*Yu20 + TYd11*Yu21 + TYd12*Yu22) + 
      TYd20*(TYd20*Yu20 + TYd21*Yu21 + TYd22*Yu22)) + 
    Yu21*(TYd01*(TYd00*Yu20 + TYd01*Yu21 + TYd02*Yu22) + 
      TYd11*(TYd10*Yu20 + TYd11*Yu21 + TYd12*Yu22) + 
      TYd21*(TYd20*Yu20 + TYd21*Yu21 + TYd22*Yu22)) + 
    Yu22*(TYd02*(TYd00*Yu20 + TYd01*Yu21 + TYd02*Yu22) + 
      TYd12*(TYd10*Yu20 + TYd11*Yu21 + TYd12*Yu22) + 
      TYd22*(TYd20*Yu20 + TYd21*Yu21 + TYd22*Yu22))) - 
  6*(TYu00*(TYd00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      TYd10*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      TYd20*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02)) + 
    TYu01*(TYd01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      TYd11*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      TYd21*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02)) + 
    TYu02*(TYd02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      TYd12*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      TYd22*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02)) + 
    TYu10*(TYd00*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      TYd10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      TYd20*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12)) + 
    TYu11*(TYd01*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      TYd11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      TYd21*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12)) + 
    TYu12*(TYd02*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      TYd12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      TYd22*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12)) + 
    TYu20*(TYd00*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 
      TYd10*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 
      TYd20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
    TYu21*(TYd01*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 
      TYd11*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 
      TYd21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
    TYu22*(TYd02*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22) + 
      TYd12*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22) + 
      TYd22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) - 
  6*(TYd00*(TYu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      TYu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      TYu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
    TYd01*(TYu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      TYu11*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      TYu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
    TYd02*(TYu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      TYu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      TYu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
    TYd10*(TYu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      TYu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      TYu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
    TYd11*(TYu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      TYu11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      TYu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
    TYd12*(TYu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      TYu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      TYu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
    TYd20*(TYu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      TYu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      TYu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
    TYd21*(TYu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      TYu11*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      TYu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
    TYd22*(TYu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      TYu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      TYu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) - 
  6*mHd2*(Yd00*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      Yu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      Yu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
    Yd01*(Yu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      Yu11*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
    Yd02*(Yu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      Yu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
    Yd10*(Yu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      Yu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      Yu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
    Yd11*(Yu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      Yu11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      Yu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
    Yd12*(Yu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      Yu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      Yu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
    Yd20*(Yu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      Yu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      Yu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
    Yd21*(Yu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      Yu11*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      Yu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
    Yd22*(Yu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      Yu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) - 
  6*mHu2*(Yd00*(Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      Yu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      Yu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
    Yd01*(Yu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      Yu11*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
    Yd02*(Yu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      Yu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
    Yd10*(Yu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      Yu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      Yu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
    Yd11*(Yu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      Yu11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      Yu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
    Yd12*(Yu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      Yu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      Yu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
    Yd20*(Yu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      Yu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      Yu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
    Yd21*(Yu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      Yu11*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      Yu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
    Yd22*(Yu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      Yu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) - 
  6*(Yd00*(Yu00*((md200*Yd00 + md201*Yd10 + md202*Yd20)*Yu00 + 
        (md200*Yd01 + md201*Yd11 + md202*Yd21)*Yu01 + 
        (md200*Yd02 + md201*Yd12 + md202*Yd22)*Yu02) + 
      Yu10*((md200*Yd00 + md201*Yd10 + md202*Yd20)*Yu10 + 
        (md200*Yd01 + md201*Yd11 + md202*Yd21)*Yu11 + 
        (md200*Yd02 + md201*Yd12 + md202*Yd22)*Yu12) + 
      Yu20*((md200*Yd00 + md201*Yd10 + md202*Yd20)*Yu20 + 
        (md200*Yd01 + md201*Yd11 + md202*Yd21)*Yu21 + 
        (md200*Yd02 + md201*Yd12 + md202*Yd22)*Yu22)) + 
    Yd01*(Yu01*((md200*Yd00 + md201*Yd10 + md202*Yd20)*Yu00 + 
        (md200*Yd01 + md201*Yd11 + md202*Yd21)*Yu01 + 
        (md200*Yd02 + md201*Yd12 + md202*Yd22)*Yu02) + 
      Yu11*((md200*Yd00 + md201*Yd10 + md202*Yd20)*Yu10 + 
        (md200*Yd01 + md201*Yd11 + md202*Yd21)*Yu11 + 
        (md200*Yd02 + md201*Yd12 + md202*Yd22)*Yu12) + 
      Yu21*((md200*Yd00 + md201*Yd10 + md202*Yd20)*Yu20 + 
        (md200*Yd01 + md201*Yd11 + md202*Yd21)*Yu21 + 
        (md200*Yd02 + md201*Yd12 + md202*Yd22)*Yu22)) + 
    Yd02*(Yu02*((md200*Yd00 + md201*Yd10 + md202*Yd20)*Yu00 + 
        (md200*Yd01 + md201*Yd11 + md202*Yd21)*Yu01 + 
        (md200*Yd02 + md201*Yd12 + md202*Yd22)*Yu02) + 
      Yu12*((md200*Yd00 + md201*Yd10 + md202*Yd20)*Yu10 + 
        (md200*Yd01 + md201*Yd11 + md202*Yd21)*Yu11 + 
        (md200*Yd02 + md201*Yd12 + md202*Yd22)*Yu12) + 
      Yu22*((md200*Yd00 + md201*Yd10 + md202*Yd20)*Yu20 + 
        (md200*Yd01 + md201*Yd11 + md202*Yd21)*Yu21 + 
        (md200*Yd02 + md201*Yd12 + md202*Yd22)*Yu22)) + 
    Yd10*(Yu00*((md210*Yd00 + md211*Yd10 + md212*Yd20)*Yu00 + 
        (md210*Yd01 + md211*Yd11 + md212*Yd21)*Yu01 + 
        (md210*Yd02 + md211*Yd12 + md212*Yd22)*Yu02) + 
      Yu10*((md210*Yd00 + md211*Yd10 + md212*Yd20)*Yu10 + 
        (md210*Yd01 + md211*Yd11 + md212*Yd21)*Yu11 + 
        (md210*Yd02 + md211*Yd12 + md212*Yd22)*Yu12) + 
      Yu20*((md210*Yd00 + md211*Yd10 + md212*Yd20)*Yu20 + 
        (md210*Yd01 + md211*Yd11 + md212*Yd21)*Yu21 + 
        (md210*Yd02 + md211*Yd12 + md212*Yd22)*Yu22)) + 
    Yd11*(Yu01*((md210*Yd00 + md211*Yd10 + md212*Yd20)*Yu00 + 
        (md210*Yd01 + md211*Yd11 + md212*Yd21)*Yu01 + 
        (md210*Yd02 + md211*Yd12 + md212*Yd22)*Yu02) + 
      Yu11*((md210*Yd00 + md211*Yd10 + md212*Yd20)*Yu10 + 
        (md210*Yd01 + md211*Yd11 + md212*Yd21)*Yu11 + 
        (md210*Yd02 + md211*Yd12 + md212*Yd22)*Yu12) + 
      Yu21*((md210*Yd00 + md211*Yd10 + md212*Yd20)*Yu20 + 
        (md210*Yd01 + md211*Yd11 + md212*Yd21)*Yu21 + 
        (md210*Yd02 + md211*Yd12 + md212*Yd22)*Yu22)) + 
    Yd12*(Yu02*((md210*Yd00 + md211*Yd10 + md212*Yd20)*Yu00 + 
        (md210*Yd01 + md211*Yd11 + md212*Yd21)*Yu01 + 
        (md210*Yd02 + md211*Yd12 + md212*Yd22)*Yu02) + 
      Yu12*((md210*Yd00 + md211*Yd10 + md212*Yd20)*Yu10 + 
        (md210*Yd01 + md211*Yd11 + md212*Yd21)*Yu11 + 
        (md210*Yd02 + md211*Yd12 + md212*Yd22)*Yu12) + 
      Yu22*((md210*Yd00 + md211*Yd10 + md212*Yd20)*Yu20 + 
        (md210*Yd01 + md211*Yd11 + md212*Yd21)*Yu21 + 
        (md210*Yd02 + md211*Yd12 + md212*Yd22)*Yu22)) + 
    Yd20*(Yu00*((md220*Yd00 + md221*Yd10 + md222*Yd20)*Yu00 + 
        (md220*Yd01 + md221*Yd11 + md222*Yd21)*Yu01 + 
        (md220*Yd02 + md221*Yd12 + md222*Yd22)*Yu02) + 
      Yu10*((md220*Yd00 + md221*Yd10 + md222*Yd20)*Yu10 + 
        (md220*Yd01 + md221*Yd11 + md222*Yd21)*Yu11 + 
        (md220*Yd02 + md221*Yd12 + md222*Yd22)*Yu12) + 
      Yu20*((md220*Yd00 + md221*Yd10 + md222*Yd20)*Yu20 + 
        (md220*Yd01 + md221*Yd11 + md222*Yd21)*Yu21 + 
        (md220*Yd02 + md221*Yd12 + md222*Yd22)*Yu22)) + 
    Yd21*(Yu01*((md220*Yd00 + md221*Yd10 + md222*Yd20)*Yu00 + 
        (md220*Yd01 + md221*Yd11 + md222*Yd21)*Yu01 + 
        (md220*Yd02 + md221*Yd12 + md222*Yd22)*Yu02) + 
      Yu11*((md220*Yd00 + md221*Yd10 + md222*Yd20)*Yu10 + 
        (md220*Yd01 + md221*Yd11 + md222*Yd21)*Yu11 + 
        (md220*Yd02 + md221*Yd12 + md222*Yd22)*Yu12) + 
      Yu21*((md220*Yd00 + md221*Yd10 + md222*Yd20)*Yu20 + 
        (md220*Yd01 + md221*Yd11 + md222*Yd21)*Yu21 + 
        (md220*Yd02 + md221*Yd12 + md222*Yd22)*Yu22)) + 
    Yd22*(Yu02*((md220*Yd00 + md221*Yd10 + md222*Yd20)*Yu00 + 
        (md220*Yd01 + md221*Yd11 + md222*Yd21)*Yu01 + 
        (md220*Yd02 + md221*Yd12 + md222*Yd22)*Yu02) + 
      Yu12*((md220*Yd00 + md221*Yd10 + md222*Yd20)*Yu10 + 
        (md220*Yd01 + md221*Yd11 + md222*Yd21)*Yu11 + 
        (md220*Yd02 + md221*Yd12 + md222*Yd22)*Yu12) + 
      Yu22*((md220*Yd00 + md221*Yd10 + md222*Yd20)*Yu20 + 
        (md220*Yd01 + md221*Yd11 + md222*Yd21)*Yu21 + 
        (md220*Yd02 + md221*Yd12 + md222*Yd22)*Yu22))) - 
  (g1^2*(-4*g1^2*(mDx200 + mDx211 + mDx222) - 
     80*g3^2*(mDx200 + mDx211 + mDx222) + 
     30*(Kappa00*(Kappa00*mDx200 + Kappa10*mDx201 + Kappa20*mDx202) + 
       Kappa01*(Kappa01*mDx200 + Kappa11*mDx201 + Kappa21*mDx202) + 
       Kappa02*(Kappa02*mDx200 + Kappa12*mDx201 + Kappa22*mDx202) + 
       Kappa10*(Kappa00*mDx210 + Kappa10*mDx211 + Kappa20*mDx212) + 
       Kappa11*(Kappa01*mDx210 + Kappa11*mDx211 + Kappa21*mDx212) + 
       Kappa12*(Kappa02*mDx210 + Kappa12*mDx211 + Kappa22*mDx212) + 
       Kappa20*(Kappa00*mDx220 + Kappa10*mDx221 + Kappa20*mDx222) + 
       Kappa21*(Kappa01*mDx220 + Kappa11*mDx221 + Kappa21*mDx222) + 
       Kappa22*(Kappa02*mDx220 + Kappa12*mDx221 + Kappa22*mDx222)) + 
     4*g1^2*(mDxbar200 + mDxbar211 + mDxbar222) + 
     80*g3^2*(mDxbar200 + mDxbar211 + mDxbar222) - 
     30*(Kappa00*(Kappa00*mDxbar200 + Kappa01*mDxbar201 + 
         Kappa02*mDxbar202) + Kappa10*(Kappa10*mDxbar200 + 
         Kappa11*mDxbar201 + Kappa12*mDxbar202) + 
       Kappa20*(Kappa20*mDxbar200 + Kappa21*mDxbar201 + Kappa22*mDxbar202) + 
       Kappa01*(Kappa00*mDxbar210 + Kappa01*mDxbar211 + Kappa02*mDxbar212) + 
       Kappa11*(Kappa10*mDxbar210 + Kappa11*mDxbar211 + Kappa12*mDxbar212) + 
       Kappa21*(Kappa20*mDxbar210 + Kappa21*mDxbar211 + Kappa22*mDxbar212) + 
       Kappa02*(Kappa00*mDxbar220 + Kappa01*mDxbar221 + Kappa02*mDxbar222) + 
       Kappa12*(Kappa10*mDxbar220 + Kappa11*mDxbar221 + Kappa12*mDxbar222) + 
       Kappa22*(Kappa20*mDxbar220 + Kappa21*mDxbar221 + Kappa22*mDxbar222)) + 
     36*g1^2*(me200 + me211 + me222) - 9*g1^2*(mH1I200 + mH1I211) - 
     45*g2^2*(mH1I200 + mH1I211) + 
     30*(Lambda1200*(Lambda1200*mH1I200 + Lambda1201*mH1I210) + 
       Lambda1210*(Lambda1210*mH1I200 + Lambda1211*mH1I210) + 
       Lambda1201*(Lambda1200*mH1I201 + Lambda1201*mH1I211) + 
       Lambda1211*(Lambda1210*mH1I201 + Lambda1211*mH1I211)) + 
     9*g1^2*(mH2I200 + mH2I211) + 45*g2^2*(mH2I200 + mH2I211) - 
     30*(Lambda1200*(Lambda1200*mH2I200 + Lambda1210*mH2I201) + 
       Lambda1201*(Lambda1201*mH2I200 + Lambda1211*mH2I201) + 
       Lambda1210*(Lambda1200*mH2I210 + Lambda1210*mH2I211) + 
       Lambda1211*(Lambda1201*mH2I210 + Lambda1211*mH2I211)) - 9*g1^2*mHd2 - 
     45*g2^2*mHd2 - 9*g1^2*mHp2 - 45*g2^2*mHp2 + 9*g1^2*mHpbar2 + 
     45*g2^2*mHpbar2 + 30*Lambdax^2*(mHd2 - mHu2) + 9*g1^2*mHu2 + 
     45*g2^2*mHu2 - 9*g1^2*(ml200 + ml211 + ml222) - 
     45*g2^2*(ml200 + ml211 + ml222) + g1^2*(mq200 + mq211 + mq222) + 
     45*g2^2*(mq200 + mq211 + mq222) + 80*g3^2*(mq200 + mq211 + mq222) - 
     32*g1^2*(mu200 + mu211 + mu222) - 160*g3^2*(mu200 + mu211 + mu222) + 
     4*(md200 + md211 + md222)*(g1^2 + 20*g3^2 + 15*gN^2*Qdp^2) + 
     60*gN^2*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp^2 - 
     60*gN^2*(mDx200 + mDx211 + mDx222)*QDxp^2 + 
     60*gN^2*(me200 + me211 + me222)*Qep^2 - 60*gN^2*(mH1I200 + mH1I211)*
      QH1p^2 - 60*gN^2*mHd2*QH1p^2 + 60*gN^2*(mH2I200 + mH2I211)*QH2p^2 + 
     60*gN^2*mHu2*QH2p^2 + 60*gN^2*mHpbar2*QHpbarp^2 - 60*gN^2*mHp2*QHpp^2 - 
     60*gN^2*(ml200 + ml211 + ml222)*QLp^2 + 60*gN^2*(mq200 + mq211 + mq222)*
      QQp^2 - 120*gN^2*(mu200 + mu211 + mu222)*Qup^2 + 
     90*mHd2*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
       Yd21^2 + Yd22^2) - 30*(Yd00*(mq200*Yd00 + mq210*Yd01 + mq220*Yd02) + 
       Yd01*(mq201*Yd00 + mq211*Yd01 + mq221*Yd02) + 
       Yd02*(mq202*Yd00 + mq212*Yd01 + mq222*Yd02) + 
       Yd10*(mq200*Yd10 + mq210*Yd11 + mq220*Yd12) + 
       Yd11*(mq201*Yd10 + mq211*Yd11 + mq221*Yd12) + 
       Yd12*(mq202*Yd10 + mq212*Yd11 + mq222*Yd12) + 
       Yd20*(mq200*Yd20 + mq210*Yd21 + mq220*Yd22) + 
       Yd21*(mq201*Yd20 + mq211*Yd21 + mq221*Yd22) + 
       Yd22*(mq202*Yd20 + mq212*Yd21 + mq222*Yd22)) - 
     60*(md200*(Yd00^2 + Yd01^2 + Yd02^2) + md201*(Yd00*Yd10 + Yd01*Yd11 + 
         Yd02*Yd12) + md210*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
       md211*(Yd10^2 + Yd11^2 + Yd12^2) + md202*(Yd00*Yd20 + Yd01*Yd21 + 
         Yd02*Yd22) + md220*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
       md212*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
       md221*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
       md222*(Yd20^2 + Yd21^2 + Yd22^2)) + 
     30*mHd2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
       Ye21^2 + Ye22^2) + 30*(Ye00*(ml200*Ye00 + ml210*Ye01 + ml220*Ye02) + 
       Ye01*(ml201*Ye00 + ml211*Ye01 + ml221*Ye02) + 
       Ye02*(ml202*Ye00 + ml212*Ye01 + ml222*Ye02) + 
       Ye10*(ml200*Ye10 + ml210*Ye11 + ml220*Ye12) + 
       Ye11*(ml201*Ye10 + ml211*Ye11 + ml221*Ye12) + 
       Ye12*(ml202*Ye10 + ml212*Ye11 + ml222*Ye12) + 
       Ye20*(ml200*Ye20 + ml210*Ye21 + ml220*Ye22) + 
       Ye21*(ml201*Ye20 + ml211*Ye21 + ml221*Ye22) + 
       Ye22*(ml202*Ye20 + ml212*Ye21 + ml222*Ye22)) - 
     60*(me200*(Ye00^2 + Ye01^2 + Ye02^2) + me201*(Ye00*Ye10 + Ye01*Ye11 + 
         Ye02*Ye12) + me210*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
       me211*(Ye10^2 + Ye11^2 + Ye12^2) + me202*(Ye00*Ye20 + Ye01*Ye21 + 
         Ye02*Ye22) + me220*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
       me212*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
       me221*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
       me222*(Ye20^2 + Ye21^2 + Ye22^2)) - 
     90*mHu2*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + Yu11^2 + Yu12^2 + Yu20^2 + 
       Yu21^2 + Yu22^2) - 30*(Yu00*(mq200*Yu00 + mq210*Yu01 + mq220*Yu02) + 
       Yu01*(mq201*Yu00 + mq211*Yu01 + mq221*Yu02) + 
       Yu02*(mq202*Yu00 + mq212*Yu01 + mq222*Yu02) + 
       Yu10*(mq200*Yu10 + mq210*Yu11 + mq220*Yu12) + 
       Yu11*(mq201*Yu10 + mq211*Yu11 + mq221*Yu12) + 
       Yu12*(mq202*Yu10 + mq212*Yu11 + mq222*Yu12) + 
       Yu20*(mq200*Yu20 + mq210*Yu21 + mq220*Yu22) + 
       Yu21*(mq201*Yu20 + mq211*Yu21 + mq221*Yu22) + 
       Yu22*(mq202*Yu20 + mq212*Yu21 + mq222*Yu22)) + 
     120*(mu200*(Yu00^2 + Yu01^2 + Yu02^2) + mu201*(Yu00*Yu10 + Yu01*Yu11 + 
         Yu02*Yu12) + mu210*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + 
       mu211*(Yu10^2 + Yu11^2 + Yu12^2) + mu202*(Yu00*Yu20 + Yu01*Yu21 + 
         Yu02*Yu22) + mu220*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22) + 
       mu212*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 
       mu221*(Yu10*Yu20 + Yu11*Yu21 + Yu12*Yu22) + 
       mu222*(Yu20^2 + Yu21^2 + Yu22^2))))/25 + 
  (4*gN^2*QH1p*(2*(md200 + md211 + md222)*Qdp*(g1^2 + 20*g3^2 + 
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
     15*g2^2*mHpbar2*QHpbarp + 20*gN^2*mHpbar2*QHpbarp^3 + 3*g1^2*mHp2*QHpp + 
     15*g2^2*mHp2*QHpp + 20*gN^2*mHp2*QHpp^3 + 3*g1^2*(ml200 + ml211 + ml222)*
      QLp + 15*g2^2*(ml200 + ml211 + ml222)*QLp + 
     20*gN^2*(ml200 + ml211 + ml222)*QLp^3 + g1^2*(mq200 + mq211 + mq222)*
      QQp + 45*g2^2*(mq200 + mq211 + mq222)*QQp + 
     80*g3^2*(mq200 + mq211 + mq222)*QQp + 60*gN^2*(mq200 + mq211 + mq222)*
      QQp^3 - 15*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
       Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*ms2*QSp - 
     10*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*ms2*QSp + 
     10*gN^2*ms2*QSp^3 + 10*gN^2*(msI200 + msI211)*QSp^3 - 
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
       mu222*(Yu20^2 + Yu21^2 + Yu22^2))))/5 - 
  6*(Yd00*(Yd00*(Yu00*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
        Yu10*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
        Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22)) + 
      Yd01*(Yu01*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
        Yu11*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
        Yu21*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22)) + 
      Yd02*(Yu02*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
        Yu12*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
        Yu22*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22))) + 
    Yd10*(Yd10*(Yu00*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
        Yu10*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
        Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22)) + 
      Yd11*(Yu01*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
        Yu11*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
        Yu21*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22)) + 
      Yd12*(Yu02*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
        Yu12*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
        Yu22*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22))) + 
    Yd20*(Yd20*(Yu00*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
        Yu10*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
        Yu20*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22)) + 
      Yd21*(Yu01*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
        Yu11*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
        Yu21*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22)) + 
      Yd22*(Yu02*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
        Yu12*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
        Yu22*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22))) + 
    Yd01*(Yd00*(Yu00*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
        Yu10*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
        Yu20*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22)) + 
      Yd01*(Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
        Yu11*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
        Yu21*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22)) + 
      Yd02*(Yu02*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
        Yu12*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
        Yu22*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22))) + 
    Yd11*(Yd10*(Yu00*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
        Yu10*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
        Yu20*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22)) + 
      Yd11*(Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
        Yu11*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
        Yu21*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22)) + 
      Yd12*(Yu02*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
        Yu12*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
        Yu22*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22))) + 
    Yd21*(Yd20*(Yu00*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
        Yu10*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
        Yu20*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22)) + 
      Yd21*(Yu01*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
        Yu11*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
        Yu21*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22)) + 
      Yd22*(Yu02*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
        Yu12*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
        Yu22*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22))) + 
    Yd02*(Yd00*(Yu00*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
        Yu10*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
        Yu20*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) + 
      Yd01*(Yu01*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
        Yu11*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
        Yu21*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) + 
      Yd02*(Yu02*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
        Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
        Yu22*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22))) + 
    Yd12*(Yd10*(Yu00*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
        Yu10*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
        Yu20*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) + 
      Yd11*(Yu01*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
        Yu11*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
        Yu21*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) + 
      Yd12*(Yu02*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
        Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
        Yu22*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22))) + 
    Yd22*(Yd20*(Yu00*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
        Yu10*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
        Yu20*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) + 
      Yd21*(Yu01*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
        Yu11*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
        Yu21*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)) + 
      Yd22*(Yu02*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
        Yu12*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
        Yu22*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22)))) - 
  6*(Yu00*(Yd00*(Yd00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
        Yd01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
        Yd02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22)) + 
      Yd10*(Yd10*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
        Yd11*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
        Yd12*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22)) + 
      Yd20*(Yd20*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
        Yd21*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
        Yd22*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22))) + 
    Yu01*(Yd01*(Yd00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
        Yd01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
        Yd02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22)) + 
      Yd11*(Yd10*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
        Yd11*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
        Yd12*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22)) + 
      Yd21*(Yd20*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
        Yd21*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
        Yd22*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22))) + 
    Yu02*(Yd02*(Yd00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
        Yd01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
        Yd02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22)) + 
      Yd12*(Yd10*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
        Yd11*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
        Yd12*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22)) + 
      Yd22*(Yd20*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
        Yd21*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
        Yd22*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22))) + 
    Yu10*(Yd00*(Yd00*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
        Yd01*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
        Yd02*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22)) + 
      Yd10*(Yd10*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
        Yd11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
        Yd12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22)) + 
      Yd20*(Yd20*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
        Yd21*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
        Yd22*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22))) + 
    Yu11*(Yd01*(Yd00*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
        Yd01*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
        Yd02*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22)) + 
      Yd11*(Yd10*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
        Yd11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
        Yd12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22)) + 
      Yd21*(Yd20*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
        Yd21*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
        Yd22*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22))) + 
    Yu12*(Yd02*(Yd00*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
        Yd01*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
        Yd02*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22)) + 
      Yd12*(Yd10*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
        Yd11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
        Yd12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22)) + 
      Yd22*(Yd20*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
        Yd21*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
        Yd22*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22))) + 
    Yu20*(Yd00*(Yd00*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
        Yd01*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
        Yd02*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) + 
      Yd10*(Yd10*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
        Yd11*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
        Yd12*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) + 
      Yd20*(Yd20*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
        Yd21*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
        Yd22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22))) + 
    Yu21*(Yd01*(Yd00*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
        Yd01*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
        Yd02*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) + 
      Yd11*(Yd10*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
        Yd11*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
        Yd12*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) + 
      Yd21*(Yd20*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
        Yd21*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
        Yd22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22))) + 
    Yu22*(Yd02*(Yd00*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
        Yd01*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
        Yd02*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) + 
      Yd12*(Yd10*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
        Yd11*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
        Yd12*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)) + 
      Yd22*(Yd20*(mu220*Yu00 + mu221*Yu10 + mu222*Yu20) + 
        Yd21*(mu220*Yu01 + mu221*Yu11 + mu222*Yu21) + 
        Yd22*(mu220*Yu02 + mu221*Yu12 + mu222*Yu22)))), 
 (-1152*g1^4*MassB^2)/25 - 144*g2^4*MassWB^2 + 
  (48*g1^3*((-12*g1*MassB^2)/5 - (6*g1*(md200 + md211 + md222 - mDx200 - 
        mDx211 - mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + 
        me222 - mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + 
        mHpbar2 + mHu2 - ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 
        2*(mu200 + mu211 + mu222)))/5))/5 - 32*gN^4*MassBp^2*QH1p^2*
   (9*Qdp^2 + 9*QDxbarp^2 + 9*QDxp^2 + 3*Qep^2 + 6*QH1p^2 + 6*QH2p^2 + 
    2*QHpbarp^2 + 2*QHpp^2 + 6*QLp^2 + 18*QQp^2 + 3*QSp^2 + 9*Qup^2) + 
  8*gN^3*QH1p*QSp^2*(-4*gN*MassBp^2*QSp + 
    gN*(3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
       QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + 
      (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 
      2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 
      2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 
      6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 
      3*(mu200 + mu211 + mu222)*Qup)) + 
  gN^3*(9*Qdp^2 + 9*QDxbarp^2 + 9*QDxp^2 + 3*Qep^2 + 6*QH1p^2 + 6*QH2p^2 + 
    2*QHpbarp^2 + 2*QHpp^2 + 6*QLp^2 + 18*QQp^2 + 3*QSp^2 + 9*Qup^2)*
   (-16*gN*MassBp^2*QH1p^2 + 4*gN*QH1p*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)) + 
  ((-3*g1^2)/5 + 4*gN^2*QH1p*QHpbarp)*((-6*g1^2*MassB^2)/5 - 
    6*g2^2*MassWB^2 + (3*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - 
       mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - 
       mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - 
       ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 
       2*(mu200 + mu211 + mu222)))/5 - 8*gN^2*MassBp^2*QHpbarp^2 + 
    2*gN^2*QHpbarp*(3*(md200 + md211 + md222)*Qdp + 
      3*(mDxbar200 + mDxbar211 + mDxbar222)*QDxbarp + 
      3*(mDx200 + mDx211 + mDx222)*QDxp + (me200 + me211 + me222)*Qep + 
      2*(mH1I200 + mH1I211)*QH1p + 2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 
      2*mHu2*QH2p + 2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 
      2*(ml200 + ml211 + ml222)*QLp + 6*(mq200 + mq211 + mq222)*QQp + 
      ms2*QSp + (msI200 + msI211)*QSp + 3*(mu200 + mu211 + mu222)*Qup)) + 
  ((3*g1^2)/5 + 4*gN^2*QH1p*QHpp)*((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 - 
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
  ((3*g1^2)/5 + 6*gN^2*QDxp*QH1p)*((-8*g1^2*MassB^2)/15 - 
    (32*g3^2*MassG^2)/3 + (Kappa00^2 + Kappa01^2 + Kappa02^2)*mDx200 + 
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
  ((3*g1^2)/5 + 6*gN^2*QDxp*QH1p)*((-8*g1^2*MassB^2)/15 - 
    (32*g3^2*MassG^2)/3 + (Kappa00*Kappa10 + Kappa01*Kappa11 + 
      Kappa02*Kappa12)*mDx201 + (Kappa10^2 + Kappa11^2 + Kappa12^2)*mDx211 + 
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
  ((-3*g1^2)/5 + 6*gN^2*QDxbarp*QH1p)*((-8*g1^2*MassB^2)/15 - 
    (32*g3^2*MassG^2)/3 + 2*(Kappa00*(Kappa00*mDx200 + Kappa10*mDx210 + 
        Kappa20*mDx220) + Kappa10*(Kappa00*mDx201 + Kappa10*mDx211 + 
        Kappa20*mDx221) + Kappa20*(Kappa00*mDx202 + Kappa10*mDx212 + 
        Kappa20*mDx222)) + (Kappa00^2 + Kappa10^2 + Kappa20^2)*mDxbar200 + 
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
  ((-3*g1^2)/5 + 6*gN^2*QDxbarp*QH1p)*((-8*g1^2*MassB^2)/15 - 
    (32*g3^2*MassG^2)/3 + 2*(Kappa01*(Kappa01*mDx200 + Kappa11*mDx210 + 
        Kappa21*mDx220) + Kappa11*(Kappa01*mDx201 + Kappa11*mDx211 + 
        Kappa21*mDx221) + Kappa21*(Kappa01*mDx202 + Kappa11*mDx212 + 
        Kappa21*mDx222)) + (Kappa00*Kappa01 + Kappa10*Kappa11 + 
      Kappa20*Kappa21)*mDxbar201 + (Kappa01^2 + Kappa11^2 + Kappa21^2)*
     mDxbar211 + Kappa01*(Kappa00*mDxbar210 + Kappa01*mDxbar211 + 
      Kappa02*mDxbar212) + Kappa11*(Kappa10*mDxbar210 + Kappa11*mDxbar211 + 
      Kappa12*mDxbar212) + Kappa21*(Kappa20*mDxbar210 + Kappa21*mDxbar211 + 
      Kappa22*mDxbar212) + (Kappa01*Kappa02 + Kappa11*Kappa12 + 
      Kappa21*Kappa22)*mDxbar221 + 2*(Kappa01^2 + Kappa11^2 + Kappa21^2)*
     ms2 + (2*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + 
       mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - 
       mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - 
       ml211 - ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/
     5 - 8*gN^2*MassBp^2*QDxbarp^2 + 2*gN^2*QDxbarp*
     (3*(md200 + md211 + md222)*Qdp + 3*(mDxbar200 + mDxbar211 + mDxbar222)*
       QDxbarp + 3*(mDx200 + mDx211 + mDx222)*QDxp + 
      (me200 + me211 + me222)*Qep + 2*(mH1I200 + mH1I211)*QH1p + 
      2*mHd2*QH1p + 2*(mH2I200 + mH2I211)*QH2p + 2*mHu2*QH2p + 
      2*mHpbar2*QHpbarp + 2*mHp2*QHpp + 2*(ml200 + ml211 + ml222)*QLp + 
      6*(mq200 + mq211 + mq222)*QQp + ms2*QSp + (msI200 + msI211)*QSp + 
      3*(mu200 + mu211 + mu222)*Qup) + 
    2*(TKappa01^2 + TKappa11^2 + TKappa21^2)) + 
  ((-3*g1^2)/5 + 6*gN^2*QDxbarp*QH1p)*((-8*g1^2*MassB^2)/15 - 
    (32*g3^2*MassG^2)/3 + 2*(Kappa02*(Kappa02*mDx200 + Kappa12*mDx210 + 
        Kappa22*mDx220) + Kappa12*(Kappa02*mDx201 + Kappa12*mDx211 + 
        Kappa22*mDx221) + Kappa22*(Kappa02*mDx202 + Kappa12*mDx212 + 
        Kappa22*mDx222)) + (Kappa00*Kappa02 + Kappa10*Kappa12 + 
      Kappa20*Kappa22)*mDxbar202 + (Kappa01*Kappa02 + Kappa11*Kappa12 + 
      Kappa21*Kappa22)*mDxbar212 + (Kappa02^2 + Kappa12^2 + Kappa22^2)*
     mDxbar222 + Kappa02*(Kappa00*mDxbar220 + Kappa01*mDxbar221 + 
      Kappa02*mDxbar222) + Kappa12*(Kappa10*mDxbar220 + Kappa11*mDxbar221 + 
      Kappa12*mDxbar222) + Kappa22*(Kappa20*mDxbar220 + Kappa21*mDxbar221 + 
      Kappa22*mDxbar222) + 2*(Kappa02^2 + Kappa12^2 + Kappa22^2)*ms2 + 
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
  ((3*g1^2)/5 + 6*gN^2*QDxp*QH1p)*((-8*g1^2*MassB^2)/15 - 
    (32*g3^2*MassG^2)/3 + (Kappa00*Kappa20 + Kappa01*Kappa21 + 
      Kappa02*Kappa22)*mDx202 + (Kappa10*Kappa20 + Kappa11*Kappa21 + 
      Kappa12*Kappa22)*mDx212 + (Kappa20^2 + Kappa21^2 + Kappa22^2)*mDx222 + 
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
  ((-3*g1^2)/5 + 4*gN^2*QH1p*QH2p)*((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 + 
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
    2*(TLambda1200^2 + TLambda1201^2)) + ((3*g1^2)/5 + 4*gN^2*QH1p^2)*
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
    2*(TLambda1200^2 + TLambda1210^2)) + ((3*g1^2)/5 + 4*gN^2*QH1p^2)*
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
    2*(TLambda1201^2 + TLambda1211^2)) + ((-3*g1^2)/5 + 4*gN^2*QH1p*QH2p)*
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
    2*(TLambda1210^2 + TLambda1211^2)) + (2*Lambdax^2 + 2*gN^2*QH1p*QSp)*
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
    4*TLambdax^2) + ((-3*g1^2)/5 + 6*gN^2*Qdp*QH1p + 
    6*(Yd00^2 + Yd01^2 + Yd02^2))*((-8*g1^2*MassB^2)/15 - 
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
  6*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12)*
   (4*(TYd00*TYd10 + TYd01*TYd11 + TYd02*TYd12) + 
    4*mHd2*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
    4*((mq200*Yd00 + mq210*Yd01 + mq220*Yd02)*Yd10 + 
      (mq201*Yd00 + mq211*Yd01 + mq221*Yd02)*Yd11 + 
      (mq202*Yd00 + mq212*Yd01 + mq222*Yd02)*Yd12) + 
    2*(Yd10*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
      Yd11*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
      Yd12*(md200*Yd02 + md201*Yd12 + md202*Yd22)) + 
    2*(md201*(Yd00^2 + Yd01^2 + Yd02^2) + md211*(Yd00*Yd10 + Yd01*Yd11 + 
        Yd02*Yd12) + md221*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22))) + 
  6*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)*
   (4*(TYd00*TYd20 + TYd01*TYd21 + TYd02*TYd22) + 
    4*mHd2*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
    4*((mq200*Yd00 + mq210*Yd01 + mq220*Yd02)*Yd20 + 
      (mq201*Yd00 + mq211*Yd01 + mq221*Yd02)*Yd21 + 
      (mq202*Yd00 + mq212*Yd01 + mq222*Yd02)*Yd22) + 
    2*(Yd20*(md200*Yd00 + md201*Yd10 + md202*Yd20) + 
      Yd21*(md200*Yd01 + md201*Yd11 + md202*Yd21) + 
      Yd22*(md200*Yd02 + md201*Yd12 + md202*Yd22)) + 
    2*(md202*(Yd00^2 + Yd01^2 + Yd02^2) + md212*(Yd00*Yd10 + Yd01*Yd11 + 
        Yd02*Yd12) + md222*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22))) + 
  6*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12)*
   (4*(TYd00*TYd10 + TYd01*TYd11 + TYd02*TYd12) + 
    4*mHd2*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
    4*(Yd00*(mq200*Yd10 + mq210*Yd11 + mq220*Yd12) + 
      Yd01*(mq201*Yd10 + mq211*Yd11 + mq221*Yd12) + 
      Yd02*(mq202*Yd10 + mq212*Yd11 + mq222*Yd12)) + 
    2*(Yd00*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
      Yd01*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
      Yd02*(md210*Yd02 + md211*Yd12 + md212*Yd22)) + 
    2*(md200*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      md210*(Yd10^2 + Yd11^2 + Yd12^2) + md220*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22))) + ((-3*g1^2)/5 + 6*gN^2*Qdp*QH1p + 
    6*(Yd10^2 + Yd11^2 + Yd12^2))*((-8*g1^2*MassB^2)/15 - 
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
    4*(TYd10^2 + TYd11^2 + TYd12^2) + 4*mHd2*(Yd10^2 + Yd11^2 + Yd12^2) + 
    4*(Yd10*(mq200*Yd10 + mq210*Yd11 + mq220*Yd12) + 
      Yd11*(mq201*Yd10 + mq211*Yd11 + mq221*Yd12) + 
      Yd12*(mq202*Yd10 + mq212*Yd11 + mq222*Yd12)) + 
    2*(Yd10*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
      Yd11*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
      Yd12*(md210*Yd02 + md211*Yd12 + md212*Yd22)) + 
    2*(md201*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      md211*(Yd10^2 + Yd11^2 + Yd12^2) + md221*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22))) + 6*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22)*
   (4*(TYd10*TYd20 + TYd11*TYd21 + TYd12*TYd22) + 
    4*mHd2*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
    4*((mq200*Yd10 + mq210*Yd11 + mq220*Yd12)*Yd20 + 
      (mq201*Yd10 + mq211*Yd11 + mq221*Yd12)*Yd21 + 
      (mq202*Yd10 + mq212*Yd11 + mq222*Yd12)*Yd22) + 
    2*(Yd20*(md210*Yd00 + md211*Yd10 + md212*Yd20) + 
      Yd21*(md210*Yd01 + md211*Yd11 + md212*Yd21) + 
      Yd22*(md210*Yd02 + md211*Yd12 + md212*Yd22)) + 
    2*(md202*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      md212*(Yd10^2 + Yd11^2 + Yd12^2) + md222*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22))) + 6*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)*
   (4*(TYd00*TYd20 + TYd01*TYd21 + TYd02*TYd22) + 
    4*mHd2*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
    2*(Yd00*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
      Yd01*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
      Yd02*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
    4*(Yd00*(mq200*Yd20 + mq210*Yd21 + mq220*Yd22) + 
      Yd01*(mq201*Yd20 + mq211*Yd21 + mq221*Yd22) + 
      Yd02*(mq202*Yd20 + mq212*Yd21 + mq222*Yd22)) + 
    2*(md200*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      md210*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      md220*(Yd20^2 + Yd21^2 + Yd22^2))) + 
  6*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22)*
   (4*(TYd10*TYd20 + TYd11*TYd21 + TYd12*TYd22) + 
    4*mHd2*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
    2*(Yd10*(md220*Yd00 + md221*Yd10 + md222*Yd20) + 
      Yd11*(md220*Yd01 + md221*Yd11 + md222*Yd21) + 
      Yd12*(md220*Yd02 + md221*Yd12 + md222*Yd22)) + 
    4*(Yd10*(mq200*Yd20 + mq210*Yd21 + mq220*Yd22) + 
      Yd11*(mq201*Yd20 + mq211*Yd21 + mq221*Yd22) + 
      Yd12*(mq202*Yd20 + mq212*Yd21 + mq222*Yd22)) + 
    2*(md201*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      md211*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      md221*(Yd20^2 + Yd21^2 + Yd22^2))) + 
  ((-3*g1^2)/5 + 6*gN^2*Qdp*QH1p + 6*(Yd20^2 + Yd21^2 + Yd22^2))*
   ((-8*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 + 
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
  ((3*g1^2)/5 + 4*gN^2*QH1p*QLp + 2*(Ye00^2 + Ye10^2 + Ye20^2))*
   ((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 - 
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
  2*(Ye00*Ye01 + Ye10*Ye11 + Ye20*Ye21)*
   (2*(TYe00*TYe01 + TYe10*TYe11 + TYe20*TYe21) + 
    Ye01*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
    Ye11*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
    ml201*(Ye00^2 + Ye10^2 + Ye20^2) + 2*mHd2*(Ye00*Ye01 + Ye10*Ye11 + 
      Ye20*Ye21) + ml211*(Ye00*Ye01 + Ye10*Ye11 + Ye20*Ye21) + 
    2*(Ye01*(me200*Ye00 + me210*Ye10 + me220*Ye20) + 
      Ye11*(me201*Ye00 + me211*Ye10 + me221*Ye20) + 
      (me202*Ye00 + me212*Ye10 + me222*Ye20)*Ye21) + 
    Ye21*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22) + 
    ml221*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22)) + 
  2*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22)*
   (2*(TYe00*TYe02 + TYe10*TYe12 + TYe20*TYe22) + 
    Ye02*(ml200*Ye00 + ml201*Ye01 + ml202*Ye02) + 
    Ye12*(ml200*Ye10 + ml201*Ye11 + ml202*Ye12) + 
    ml202*(Ye00^2 + Ye10^2 + Ye20^2) + ml212*(Ye00*Ye01 + Ye10*Ye11 + 
      Ye20*Ye21) + Ye22*(ml200*Ye20 + ml201*Ye21 + ml202*Ye22) + 
    2*mHd2*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22) + 
    ml222*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22) + 
    2*(Ye02*(me200*Ye00 + me210*Ye10 + me220*Ye20) + 
      Ye12*(me201*Ye00 + me211*Ye10 + me221*Ye20) + 
      (me202*Ye00 + me212*Ye10 + me222*Ye20)*Ye22)) + 
  2*(Ye00*Ye01 + Ye10*Ye11 + Ye20*Ye21)*
   (2*(TYe00*TYe01 + TYe10*TYe11 + TYe20*TYe21) + 
    Ye00*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
    Ye10*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
    2*mHd2*(Ye00*Ye01 + Ye10*Ye11 + Ye20*Ye21) + 
    ml200*(Ye00*Ye01 + Ye10*Ye11 + Ye20*Ye21) + 
    ml210*(Ye01^2 + Ye11^2 + Ye21^2) + 
    2*(Ye00*(me200*Ye01 + me210*Ye11 + me220*Ye21) + 
      Ye10*(me201*Ye01 + me211*Ye11 + me221*Ye21) + 
      Ye20*(me202*Ye01 + me212*Ye11 + me222*Ye21)) + 
    Ye20*(ml210*Ye20 + ml211*Ye21 + ml212*Ye22) + 
    ml220*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22)) + 
  ((3*g1^2)/5 + 4*gN^2*QH1p*QLp + 2*(Ye01^2 + Ye11^2 + Ye21^2))*
   ((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 - 
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
  2*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22)*
   (2*(TYe01*TYe02 + TYe11*TYe12 + TYe21*TYe22) + 
    Ye02*(ml210*Ye00 + ml211*Ye01 + ml212*Ye02) + 
    Ye12*(ml210*Ye10 + ml211*Ye11 + ml212*Ye12) + 
    ml202*(Ye00*Ye01 + Ye10*Ye11 + Ye20*Ye21) + 
    ml212*(Ye01^2 + Ye11^2 + Ye21^2) + Ye22*(ml210*Ye20 + ml211*Ye21 + 
      ml212*Ye22) + 2*mHd2*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22) + 
    ml222*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22) + 
    2*(Ye02*(me200*Ye01 + me210*Ye11 + me220*Ye21) + 
      Ye12*(me201*Ye01 + me211*Ye11 + me221*Ye21) + 
      (me202*Ye01 + me212*Ye11 + me222*Ye21)*Ye22)) + 
  2*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22)*
   (2*(TYe00*TYe02 + TYe10*TYe12 + TYe20*TYe22) + 
    Ye00*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
    Ye10*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
    Ye20*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22) + 
    2*mHd2*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22) + 
    ml200*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22) + 
    ml210*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22) + 
    ml220*(Ye02^2 + Ye12^2 + Ye22^2) + 
    2*(Ye00*(me200*Ye02 + me210*Ye12 + me220*Ye22) + 
      Ye10*(me201*Ye02 + me211*Ye12 + me221*Ye22) + 
      Ye20*(me202*Ye02 + me212*Ye12 + me222*Ye22))) + 
  2*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22)*
   (2*(TYe01*TYe02 + TYe11*TYe12 + TYe21*TYe22) + 
    Ye01*(ml220*Ye00 + ml221*Ye01 + ml222*Ye02) + 
    Ye11*(ml220*Ye10 + ml221*Ye11 + ml222*Ye12) + 
    Ye21*(ml220*Ye20 + ml221*Ye21 + ml222*Ye22) + 
    ml201*(Ye00*Ye02 + Ye10*Ye12 + Ye20*Ye22) + 
    2*mHd2*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22) + 
    ml211*(Ye01*Ye02 + Ye11*Ye12 + Ye21*Ye22) + 
    ml221*(Ye02^2 + Ye12^2 + Ye22^2) + 
    2*(Ye01*(me200*Ye02 + me210*Ye12 + me220*Ye22) + 
      Ye11*(me201*Ye02 + me211*Ye12 + me221*Ye22) + 
      Ye21*(me202*Ye02 + me212*Ye12 + me222*Ye22))) + 
  ((3*g1^2)/5 + 4*gN^2*QH1p*QLp + 2*(Ye02^2 + Ye12^2 + Ye22^2))*
   ((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 - 
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
  ((3*g1^2)/5 + 2*Lambdax^2 + 4*gN^2*QH1p^2 + 
    6*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
      Yd21^2 + Yd22^2) + 2*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2))*((-6*g1^2*MassB^2)/5 - 
    6*g2^2*MassWB^2 + 2*Lambdax^2*mHd2 + 2*Lambdax^2*mHu2 + 2*Lambdax^2*ms2 - 
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
  ((-3*g1^2)/5 + 2*gN^2*Qep*QH1p + 2*(Ye00^2 + Ye01^2 + Ye02^2))*
   ((-24*g1^2*MassB^2)/5 + (6*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - 
       mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - 
       mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - 
       ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 
       2*(mu200 + mu211 + mu222)))/5 - 8*gN^2*MassBp^2*Qep^2 + 
    2*gN^2*Qep*(3*(md200 + md211 + md222)*Qdp + 
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
  2*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12)*
   (4*(TYe00*TYe10 + TYe01*TYe11 + TYe02*TYe12) + 
    4*mHd2*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
    4*((ml200*Ye00 + ml210*Ye01 + ml220*Ye02)*Ye10 + 
      (ml201*Ye00 + ml211*Ye01 + ml221*Ye02)*Ye11 + 
      (ml202*Ye00 + ml212*Ye01 + ml222*Ye02)*Ye12) + 
    2*(Ye10*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
      Ye11*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
      Ye12*(me200*Ye02 + me201*Ye12 + me202*Ye22)) + 
    2*(me201*(Ye00^2 + Ye01^2 + Ye02^2) + me211*(Ye00*Ye10 + Ye01*Ye11 + 
        Ye02*Ye12) + me221*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22))) + 
  2*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22)*
   (4*(TYe00*TYe20 + TYe01*TYe21 + TYe02*TYe22) + 
    4*mHd2*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
    4*((ml200*Ye00 + ml210*Ye01 + ml220*Ye02)*Ye20 + 
      (ml201*Ye00 + ml211*Ye01 + ml221*Ye02)*Ye21 + 
      (ml202*Ye00 + ml212*Ye01 + ml222*Ye02)*Ye22) + 
    2*(Ye20*(me200*Ye00 + me201*Ye10 + me202*Ye20) + 
      Ye21*(me200*Ye01 + me201*Ye11 + me202*Ye21) + 
      Ye22*(me200*Ye02 + me201*Ye12 + me202*Ye22)) + 
    2*(me202*(Ye00^2 + Ye01^2 + Ye02^2) + me212*(Ye00*Ye10 + Ye01*Ye11 + 
        Ye02*Ye12) + me222*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22))) + 
  4*TYe00*((-9*g1^2*TYe00)/5 - 3*g2^2*TYe00 + Lambdax^2*TYe00 - 
    2*gN^2*Qep^2*TYe00 - 2*gN^2*QH1p^2*TYe00 - 2*gN^2*QLp^2*TYe00 + 
    3*TYe00*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
      Yd21^2 + Yd22^2) + TYe00*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    5*(Ye00*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02) + 
      Ye10*(TYe00*Ye10 + TYe01*Ye11 + TYe02*Ye12) + 
      Ye20*(TYe00*Ye20 + TYe01*Ye21 + TYe02*Ye22)) + 
    Ye00*((18*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*Qep^2 + 
      4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QLp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    4*(TYe00*(Ye00^2 + Ye01^2 + Ye02^2) + TYe10*(Ye00*Ye10 + Ye01*Ye11 + 
        Ye02*Ye12) + TYe20*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22))) + 
  4*TYe01*((-9*g1^2*TYe01)/5 - 3*g2^2*TYe01 + Lambdax^2*TYe01 - 
    2*gN^2*Qep^2*TYe01 - 2*gN^2*QH1p^2*TYe01 - 2*gN^2*QLp^2*TYe01 + 
    3*TYe01*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
      Yd21^2 + Yd22^2) + TYe01*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    5*(Ye01*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02) + 
      Ye11*(TYe00*Ye10 + TYe01*Ye11 + TYe02*Ye12) + 
      Ye21*(TYe00*Ye20 + TYe01*Ye21 + TYe02*Ye22)) + 
    Ye01*((18*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*Qep^2 + 
      4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QLp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    4*(TYe01*(Ye00^2 + Ye01^2 + Ye02^2) + TYe11*(Ye00*Ye10 + Ye01*Ye11 + 
        Ye02*Ye12) + TYe21*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22))) + 
  4*TYe02*((-9*g1^2*TYe02)/5 - 3*g2^2*TYe02 + Lambdax^2*TYe02 - 
    2*gN^2*Qep^2*TYe02 - 2*gN^2*QH1p^2*TYe02 - 2*gN^2*QLp^2*TYe02 + 
    3*TYe02*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
      Yd21^2 + Yd22^2) + TYe02*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    5*(Ye02*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02) + 
      Ye12*(TYe00*Ye10 + TYe01*Ye11 + TYe02*Ye12) + 
      Ye22*(TYe00*Ye20 + TYe01*Ye21 + TYe02*Ye22)) + 
    Ye02*((18*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*Qep^2 + 
      4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QLp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    4*(TYe02*(Ye00^2 + Ye01^2 + Ye02^2) + TYe12*(Ye00*Ye10 + Ye01*Ye11 + 
        Ye02*Ye12) + TYe22*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22))) + 
  (4*mHd2*Ye00 + 2*(2*ml200*Ye00 + ml201*Ye01 + ml210*Ye01 + ml202*Ye02 + 
      ml220*Ye02) + 2*(2*me200*Ye00 + me201*Ye10 + me210*Ye10 + me202*Ye20 + 
      me220*Ye20))*(Ye00*((-9*g1^2)/5 - 3*g2^2 + Lambdax^2 - 2*gN^2*Qep^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QLp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    3*(Ye00*(Ye00^2 + Ye01^2 + Ye02^2) + Ye10*(Ye00*Ye10 + Ye01*Ye11 + 
        Ye02*Ye12) + Ye20*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22))) + 
  (4*mHd2*Ye01 + 2*(ml201*Ye00 + ml210*Ye00 + 2*ml211*Ye01 + ml212*Ye02 + 
      ml221*Ye02) + 2*(2*me200*Ye01 + me201*Ye11 + me210*Ye11 + me202*Ye21 + 
      me220*Ye21))*(Ye01*((-9*g1^2)/5 - 3*g2^2 + Lambdax^2 - 2*gN^2*Qep^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QLp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    3*(Ye01*(Ye00^2 + Ye01^2 + Ye02^2) + Ye11*(Ye00*Ye10 + Ye01*Ye11 + 
        Ye02*Ye12) + Ye21*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22))) + 
  (4*mHd2*Ye02 + 2*(ml202*Ye00 + ml220*Ye00 + ml212*Ye01 + ml221*Ye01 + 
      2*ml222*Ye02) + 2*(2*me200*Ye02 + me201*Ye12 + me210*Ye12 + 
      me202*Ye22 + me220*Ye22))*
   (Ye02*((-9*g1^2)/5 - 3*g2^2 + Lambdax^2 - 2*gN^2*Qep^2 - 2*gN^2*QH1p^2 - 
      2*gN^2*QLp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
        Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + 
      Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    3*(Ye02*(Ye00^2 + Ye01^2 + Ye02^2) + Ye12*(Ye00*Ye10 + Ye01*Ye11 + 
        Ye02*Ye12) + Ye22*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22))) + 
  2*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12)*
   (4*(TYe00*TYe10 + TYe01*TYe11 + TYe02*TYe12) + 
    4*mHd2*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
    4*(Ye00*(ml200*Ye10 + ml210*Ye11 + ml220*Ye12) + 
      Ye01*(ml201*Ye10 + ml211*Ye11 + ml221*Ye12) + 
      Ye02*(ml202*Ye10 + ml212*Ye11 + ml222*Ye12)) + 
    2*(Ye00*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
      Ye01*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
      Ye02*(me210*Ye02 + me211*Ye12 + me212*Ye22)) + 
    2*(me200*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      me210*(Ye10^2 + Ye11^2 + Ye12^2) + me220*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22))) + ((-3*g1^2)/5 + 2*gN^2*Qep*QH1p + 
    2*(Ye10^2 + Ye11^2 + Ye12^2))*((-24*g1^2*MassB^2)/5 + 
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
        Ye12*Ye22))) + 2*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22)*
   (4*(TYe10*TYe20 + TYe11*TYe21 + TYe12*TYe22) + 
    4*mHd2*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
    4*((ml200*Ye10 + ml210*Ye11 + ml220*Ye12)*Ye20 + 
      (ml201*Ye10 + ml211*Ye11 + ml221*Ye12)*Ye21 + 
      (ml202*Ye10 + ml212*Ye11 + ml222*Ye12)*Ye22) + 
    2*(Ye20*(me210*Ye00 + me211*Ye10 + me212*Ye20) + 
      Ye21*(me210*Ye01 + me211*Ye11 + me212*Ye21) + 
      Ye22*(me210*Ye02 + me211*Ye12 + me212*Ye22)) + 
    2*(me202*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      me212*(Ye10^2 + Ye11^2 + Ye12^2) + me222*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22))) + 4*TYe10*((-9*g1^2*TYe10)/5 - 3*g2^2*TYe10 + 
    Lambdax^2*TYe10 - 2*gN^2*Qep^2*TYe10 - 2*gN^2*QH1p^2*TYe10 - 
    2*gN^2*QLp^2*TYe10 + 3*TYe10*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    TYe10*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + 5*(Ye00*(TYe10*Ye00 + TYe11*Ye01 + TYe12*Ye02) + 
      Ye10*(TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12) + 
      Ye20*(TYe10*Ye20 + TYe11*Ye21 + TYe12*Ye22)) + 
    Ye10*((18*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*Qep^2 + 
      4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QLp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    4*(TYe00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      TYe10*(Ye10^2 + Ye11^2 + Ye12^2) + TYe20*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22))) + 4*TYe11*((-9*g1^2*TYe11)/5 - 3*g2^2*TYe11 + 
    Lambdax^2*TYe11 - 2*gN^2*Qep^2*TYe11 - 2*gN^2*QH1p^2*TYe11 - 
    2*gN^2*QLp^2*TYe11 + 3*TYe11*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    TYe11*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + 5*(Ye01*(TYe10*Ye00 + TYe11*Ye01 + TYe12*Ye02) + 
      Ye11*(TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12) + 
      Ye21*(TYe10*Ye20 + TYe11*Ye21 + TYe12*Ye22)) + 
    Ye11*((18*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*Qep^2 + 
      4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QLp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    4*(TYe01*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      TYe11*(Ye10^2 + Ye11^2 + Ye12^2) + TYe21*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22))) + 4*TYe12*((-9*g1^2*TYe12)/5 - 3*g2^2*TYe12 + 
    Lambdax^2*TYe12 - 2*gN^2*Qep^2*TYe12 - 2*gN^2*QH1p^2*TYe12 - 
    2*gN^2*QLp^2*TYe12 + 3*TYe12*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    TYe12*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + 5*(Ye02*(TYe10*Ye00 + TYe11*Ye01 + TYe12*Ye02) + 
      Ye12*(TYe10*Ye10 + TYe11*Ye11 + TYe12*Ye12) + 
      Ye22*(TYe10*Ye20 + TYe11*Ye21 + TYe12*Ye22)) + 
    Ye12*((18*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*Qep^2 + 
      4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QLp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    4*(TYe02*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      TYe12*(Ye10^2 + Ye11^2 + Ye12^2) + TYe22*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22))) + 
  (4*mHd2*Ye10 + 2*(2*ml200*Ye10 + ml201*Ye11 + ml210*Ye11 + ml202*Ye12 + 
      ml220*Ye12) + 2*(me201*Ye00 + me210*Ye00 + 2*me211*Ye10 + me212*Ye20 + 
      me221*Ye20))*(Ye10*((-9*g1^2)/5 - 3*g2^2 + Lambdax^2 - 2*gN^2*Qep^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QLp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    3*(Ye00*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      Ye10*(Ye10^2 + Ye11^2 + Ye12^2) + Ye20*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22))) + 
  (4*mHd2*Ye11 + 2*(ml201*Ye10 + ml210*Ye10 + 2*ml211*Ye11 + ml212*Ye12 + 
      ml221*Ye12) + 2*(me201*Ye01 + me210*Ye01 + 2*me211*Ye11 + me212*Ye21 + 
      me221*Ye21))*(Ye11*((-9*g1^2)/5 - 3*g2^2 + Lambdax^2 - 2*gN^2*Qep^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QLp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    3*(Ye01*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      Ye11*(Ye10^2 + Ye11^2 + Ye12^2) + Ye21*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22))) + 
  (4*mHd2*Ye12 + 2*(ml202*Ye10 + ml220*Ye10 + ml212*Ye11 + ml221*Ye11 + 
      2*ml222*Ye12) + 2*(me201*Ye02 + me210*Ye02 + 2*me211*Ye12 + 
      me212*Ye22 + me221*Ye22))*
   (Ye12*((-9*g1^2)/5 - 3*g2^2 + Lambdax^2 - 2*gN^2*Qep^2 - 2*gN^2*QH1p^2 - 
      2*gN^2*QLp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
        Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + 
      Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    3*(Ye02*(Ye00*Ye10 + Ye01*Ye11 + Ye02*Ye12) + 
      Ye12*(Ye10^2 + Ye11^2 + Ye12^2) + Ye22*(Ye10*Ye20 + Ye11*Ye21 + 
        Ye12*Ye22))) + 2*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22)*
   (4*(TYe00*TYe20 + TYe01*TYe21 + TYe02*TYe22) + 
    4*mHd2*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
    2*(Ye00*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
      Ye01*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
      Ye02*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
    4*(Ye00*(ml200*Ye20 + ml210*Ye21 + ml220*Ye22) + 
      Ye01*(ml201*Ye20 + ml211*Ye21 + ml221*Ye22) + 
      Ye02*(ml202*Ye20 + ml212*Ye21 + ml222*Ye22)) + 
    2*(me200*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      me210*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      me220*(Ye20^2 + Ye21^2 + Ye22^2))) + 
  2*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22)*
   (4*(TYe10*TYe20 + TYe11*TYe21 + TYe12*TYe22) + 
    4*mHd2*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
    2*(Ye10*(me220*Ye00 + me221*Ye10 + me222*Ye20) + 
      Ye11*(me220*Ye01 + me221*Ye11 + me222*Ye21) + 
      Ye12*(me220*Ye02 + me221*Ye12 + me222*Ye22)) + 
    4*(Ye10*(ml200*Ye20 + ml210*Ye21 + ml220*Ye22) + 
      Ye11*(ml201*Ye20 + ml211*Ye21 + ml221*Ye22) + 
      Ye12*(ml202*Ye20 + ml212*Ye21 + ml222*Ye22)) + 
    2*(me201*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      me211*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      me221*(Ye20^2 + Ye21^2 + Ye22^2))) + 
  ((-3*g1^2)/5 + 2*gN^2*Qep*QH1p + 2*(Ye20^2 + Ye21^2 + Ye22^2))*
   ((-24*g1^2*MassB^2)/5 + (6*g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - 
       mDx222 + mDxbar200 + mDxbar211 + mDxbar222 + me200 + me211 + me222 - 
       mH1I200 - mH1I211 + mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - 
       ml200 - ml211 - ml222 + mq200 + mq211 + mq222 - 
       2*(mu200 + mu211 + mu222)))/5 - 8*gN^2*MassBp^2*Qep^2 + 
    2*gN^2*Qep*(3*(md200 + md211 + md222)*Qdp + 
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
  4*TYe20*((-9*g1^2*TYe20)/5 - 3*g2^2*TYe20 + Lambdax^2*TYe20 - 
    2*gN^2*Qep^2*TYe20 - 2*gN^2*QH1p^2*TYe20 - 2*gN^2*QLp^2*TYe20 + 
    3*TYe20*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
      Yd21^2 + Yd22^2) + TYe20*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    5*(Ye00*(TYe20*Ye00 + TYe21*Ye01 + TYe22*Ye02) + 
      Ye10*(TYe20*Ye10 + TYe21*Ye11 + TYe22*Ye12) + 
      Ye20*(TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Ye20*((18*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*Qep^2 + 
      4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QLp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    4*(TYe00*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      TYe10*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      TYe20*(Ye20^2 + Ye21^2 + Ye22^2))) + 
  4*TYe21*((-9*g1^2*TYe21)/5 - 3*g2^2*TYe21 + Lambdax^2*TYe21 - 
    2*gN^2*Qep^2*TYe21 - 2*gN^2*QH1p^2*TYe21 - 2*gN^2*QLp^2*TYe21 + 
    3*TYe21*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
      Yd21^2 + Yd22^2) + TYe21*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    5*(Ye01*(TYe20*Ye00 + TYe21*Ye01 + TYe22*Ye02) + 
      Ye11*(TYe20*Ye10 + TYe21*Ye11 + TYe22*Ye12) + 
      Ye21*(TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Ye21*((18*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*Qep^2 + 
      4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QLp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    4*(TYe01*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      TYe11*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      TYe21*(Ye20^2 + Ye21^2 + Ye22^2))) + 
  4*TYe22*((-9*g1^2*TYe22)/5 - 3*g2^2*TYe22 + Lambdax^2*TYe22 - 
    2*gN^2*Qep^2*TYe22 - 2*gN^2*QH1p^2*TYe22 - 2*gN^2*QLp^2*TYe22 + 
    3*TYe22*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
      Yd21^2 + Yd22^2) + TYe22*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    5*(Ye02*(TYe20*Ye00 + TYe21*Ye01 + TYe22*Ye02) + 
      Ye12*(TYe20*Ye10 + TYe21*Ye11 + TYe22*Ye12) + 
      Ye22*(TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Ye22*((18*g1^2*MassB)/5 + 6*g2^2*MassWB + 4*gN^2*MassBp*Qep^2 + 
      4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QLp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    4*(TYe02*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      TYe12*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      TYe22*(Ye20^2 + Ye21^2 + Ye22^2))) + 
  (4*mHd2*Ye20 + 2*(me202*Ye00 + me220*Ye00 + me212*Ye10 + me221*Ye10 + 
      2*me222*Ye20) + 2*(2*ml200*Ye20 + ml201*Ye21 + ml210*Ye21 + 
      ml202*Ye22 + ml220*Ye22))*
   (Ye20*((-9*g1^2)/5 - 3*g2^2 + Lambdax^2 - 2*gN^2*Qep^2 - 2*gN^2*QH1p^2 - 
      2*gN^2*QLp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
        Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + 
      Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    3*(Ye00*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      Ye10*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      Ye20*(Ye20^2 + Ye21^2 + Ye22^2))) + 
  (4*mHd2*Ye21 + 2*(me202*Ye01 + me220*Ye01 + me212*Ye11 + me221*Ye11 + 
      2*me222*Ye21) + 2*(ml201*Ye20 + ml210*Ye20 + 2*ml211*Ye21 + 
      ml212*Ye22 + ml221*Ye22))*
   (Ye21*((-9*g1^2)/5 - 3*g2^2 + Lambdax^2 - 2*gN^2*Qep^2 - 2*gN^2*QH1p^2 - 
      2*gN^2*QLp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + 
        Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + 
      Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    3*(Ye01*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      Ye11*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      Ye21*(Ye20^2 + Ye21^2 + Ye22^2))) + 
  (4*mHd2*Ye22 + 2*(me202*Ye02 + me220*Ye02 + me212*Ye12 + me221*Ye12 + 
      2*me222*Ye22) + 2*(ml202*Ye20 + ml220*Ye20 + ml212*Ye21 + ml221*Ye21 + 
      2*ml222*Ye22))*(Ye22*((-9*g1^2)/5 - 3*g2^2 + Lambdax^2 - 2*gN^2*Qep^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QLp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    3*(Ye02*(Ye00*Ye20 + Ye01*Ye21 + Ye02*Ye22) + 
      Ye12*(Ye10*Ye20 + Ye11*Ye21 + Ye12*Ye22) + 
      Ye22*(Ye20^2 + Ye21^2 + Ye22^2))) + 
  (12*mHd2*Yd00 + 6*(2*mq200*Yd00 + mq201*Yd01 + mq210*Yd01 + mq202*Yd02 + 
      mq220*Yd02) + 6*(2*md200*Yd00 + md201*Yd10 + md210*Yd10 + md202*Yd20 + 
      md220*Yd20))*(3*(Yd00*(Yd00^2 + Yd01^2 + Yd02^2) + 
      Yd10*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd20*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    Yd00*((-7*g1^2)/15 - 3*g2^2 - (16*g3^2)/3 + Lambdax^2 - 2*gN^2*Qdp^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QQp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
    Yu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
    Yu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
  (12*mHd2*Yd01 + 6*(mq201*Yd00 + mq210*Yd00 + 2*mq211*Yd01 + mq212*Yd02 + 
      mq221*Yd02) + 6*(2*md200*Yd01 + md201*Yd11 + md210*Yd11 + md202*Yd21 + 
      md220*Yd21))*(3*(Yd01*(Yd00^2 + Yd01^2 + Yd02^2) + 
      Yd11*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    Yd01*((-7*g1^2)/15 - 3*g2^2 - (16*g3^2)/3 + Lambdax^2 - 2*gN^2*Qdp^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QQp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
    Yu11*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
    Yu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
  (12*mHd2*Yd02 + 6*(mq202*Yd00 + mq220*Yd00 + mq212*Yd01 + mq221*Yd01 + 
      2*mq222*Yd02) + 6*(2*md200*Yd02 + md201*Yd12 + md210*Yd12 + 
      md202*Yd22 + md220*Yd22))*
   (3*(Yd02*(Yd00^2 + Yd01^2 + Yd02^2) + Yd12*(Yd00*Yd10 + Yd01*Yd11 + 
        Yd02*Yd12) + Yd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    Yd02*((-7*g1^2)/15 - 3*g2^2 - (16*g3^2)/3 + Lambdax^2 - 2*gN^2*Qdp^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QQp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
    Yu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
    Yu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22)) + 
  (12*mHd2*Yd10 + 6*(2*mq200*Yd10 + mq201*Yd11 + mq210*Yd11 + mq202*Yd12 + 
      mq220*Yd12) + 6*(md201*Yd00 + md210*Yd00 + 2*md211*Yd10 + md212*Yd20 + 
      md221*Yd20))*(3*(Yd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd10*(Yd10^2 + Yd11^2 + Yd12^2) + Yd20*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + Yd10*((-7*g1^2)/15 - 3*g2^2 - (16*g3^2)/3 + Lambdax^2 - 
      2*gN^2*Qdp^2 - 2*gN^2*QH1p^2 - 2*gN^2*QQp^2 + 
      3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
        Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
    Yu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
    Yu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
  (12*mHd2*Yd11 + 6*(mq201*Yd10 + mq210*Yd10 + 2*mq211*Yd11 + mq212*Yd12 + 
      mq221*Yd12) + 6*(md201*Yd01 + md210*Yd01 + 2*md211*Yd11 + md212*Yd21 + 
      md221*Yd21))*(3*(Yd01*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd11*(Yd10^2 + Yd11^2 + Yd12^2) + Yd21*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + Yd11*((-7*g1^2)/15 - 3*g2^2 - (16*g3^2)/3 + Lambdax^2 - 
      2*gN^2*Qdp^2 - 2*gN^2*QH1p^2 - 2*gN^2*QQp^2 + 
      3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
        Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
    Yu11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
    Yu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
  (12*mHd2*Yd12 + 6*(mq202*Yd10 + mq220*Yd10 + mq212*Yd11 + mq221*Yd11 + 
      2*mq222*Yd12) + 6*(md201*Yd02 + md210*Yd02 + 2*md211*Yd12 + 
      md212*Yd22 + md221*Yd22))*
   (3*(Yd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      Yd12*(Yd10^2 + Yd11^2 + Yd12^2) + Yd22*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + Yd12*((-7*g1^2)/15 - 3*g2^2 - (16*g3^2)/3 + Lambdax^2 - 
      2*gN^2*Qdp^2 - 2*gN^2*QH1p^2 - 2*gN^2*QQp^2 + 
      3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + Yd11^2 + Yd12^2 + Yd20^2 + 
        Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
    Yu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
    Yu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22)) + 
  (12*mHd2*Yd20 + 6*(md202*Yd00 + md220*Yd00 + md212*Yd10 + md221*Yd10 + 
      2*md222*Yd20) + 6*(2*mq200*Yd20 + mq201*Yd21 + mq210*Yd21 + 
      mq202*Yd22 + mq220*Yd22))*
   (3*(Yd00*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      Yd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      Yd20*(Yd20^2 + Yd21^2 + Yd22^2)) + 
    Yd20*((-7*g1^2)/15 - 3*g2^2 - (16*g3^2)/3 + Lambdax^2 - 2*gN^2*Qdp^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QQp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
    Yu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
    Yu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
  (12*mHd2*Yd21 + 6*(md202*Yd01 + md220*Yd01 + md212*Yd11 + md221*Yd11 + 
      2*md222*Yd21) + 6*(mq201*Yd20 + mq210*Yd20 + 2*mq211*Yd21 + 
      mq212*Yd22 + mq221*Yd22))*
   (3*(Yd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      Yd11*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      Yd21*(Yd20^2 + Yd21^2 + Yd22^2)) + 
    Yd21*((-7*g1^2)/15 - 3*g2^2 - (16*g3^2)/3 + Lambdax^2 - 2*gN^2*Qdp^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QQp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
    Yu11*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
    Yu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
  (12*mHd2*Yd22 + 6*(md202*Yd02 + md220*Yd02 + md212*Yd12 + md221*Yd12 + 
      2*md222*Yd22) + 6*(mq202*Yd20 + mq220*Yd20 + mq212*Yd21 + mq221*Yd21 + 
      2*mq222*Yd22))*(3*(Yd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      Yd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      Yd22*(Yd20^2 + Yd21^2 + Yd22^2)) + 
    Yd22*((-7*g1^2)/15 - 3*g2^2 - (16*g3^2)/3 + Lambdax^2 - 2*gN^2*Qdp^2 - 
      2*gN^2*QH1p^2 - 2*gN^2*QQp^2 + 3*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
        Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + Ye00^2 + Ye01^2 + 
      Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
    Yu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
    Yu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22)) + 
  ((-3*g1^2)/5 + 12*gN^2*QH1p*QQp + 6*(Yd00^2 + Yd10^2 + Yd20^2))*
   ((-2*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 - 6*g2^2*MassWB^2 + 
    (g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QQp^2 + 2*gN^2*QQp*(3*(md200 + md211 + md222)*Qdp + 
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
  6*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21)*
   (2*(TYd00*TYd01 + TYd10*TYd11 + TYd20*TYd21) + 
    2*(TYu00*TYu01 + TYu10*TYu11 + TYu20*TYu21) + 
    Yd01*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
    Yd11*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
    mq201*(Yd00^2 + Yd10^2 + Yd20^2) + 2*mHd2*(Yd00*Yd01 + Yd10*Yd11 + 
      Yd20*Yd21) + mq211*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21) + 
    2*(Yd01*(md200*Yd00 + md210*Yd10 + md220*Yd20) + 
      Yd11*(md201*Yd00 + md211*Yd10 + md221*Yd20) + 
      (md202*Yd00 + md212*Yd10 + md222*Yd20)*Yd21) + 
    Yd21*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + 
    mq221*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + 
    Yu01*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
    Yu11*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
    mq201*(Yu00^2 + Yu10^2 + Yu20^2) + 2*mHu2*(Yu00*Yu01 + Yu10*Yu11 + 
      Yu20*Yu21) + mq211*(Yu00*Yu01 + Yu10*Yu11 + Yu20*Yu21) + 
    2*(Yu01*(mu200*Yu00 + mu210*Yu10 + mu220*Yu20) + 
      Yu11*(mu201*Yu00 + mu211*Yu10 + mu221*Yu20) + 
      (mu202*Yu00 + mu212*Yu10 + mu222*Yu20)*Yu21) + 
    Yu21*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22) + 
    mq221*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22)) + 
  6*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22)*
   (2*(TYd00*TYd02 + TYd10*TYd12 + TYd20*TYd22) + 
    2*(TYu00*TYu02 + TYu10*TYu12 + TYu20*TYu22) + 
    Yd02*(mq200*Yd00 + mq201*Yd01 + mq202*Yd02) + 
    Yd12*(mq200*Yd10 + mq201*Yd11 + mq202*Yd12) + 
    mq202*(Yd00^2 + Yd10^2 + Yd20^2) + mq212*(Yd00*Yd01 + Yd10*Yd11 + 
      Yd20*Yd21) + Yd22*(mq200*Yd20 + mq201*Yd21 + mq202*Yd22) + 
    2*mHd2*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + 
    mq222*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + 
    2*(Yd02*(md200*Yd00 + md210*Yd10 + md220*Yd20) + 
      Yd12*(md201*Yd00 + md211*Yd10 + md221*Yd20) + 
      (md202*Yd00 + md212*Yd10 + md222*Yd20)*Yd22) + 
    Yu02*(mq200*Yu00 + mq201*Yu01 + mq202*Yu02) + 
    Yu12*(mq200*Yu10 + mq201*Yu11 + mq202*Yu12) + 
    mq202*(Yu00^2 + Yu10^2 + Yu20^2) + mq212*(Yu00*Yu01 + Yu10*Yu11 + 
      Yu20*Yu21) + Yu22*(mq200*Yu20 + mq201*Yu21 + mq202*Yu22) + 
    2*mHu2*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 
    mq222*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 
    2*(Yu02*(mu200*Yu00 + mu210*Yu10 + mu220*Yu20) + 
      Yu12*(mu201*Yu00 + mu211*Yu10 + mu221*Yu20) + 
      (mu202*Yu00 + mu212*Yu10 + mu222*Yu20)*Yu22)) + 
  6*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21)*
   (2*(TYd00*TYd01 + TYd10*TYd11 + TYd20*TYd21) + 
    2*(TYu00*TYu01 + TYu10*TYu11 + TYu20*TYu21) + 
    Yd00*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
    Yd10*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
    2*mHd2*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21) + 
    mq200*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21) + 
    mq210*(Yd01^2 + Yd11^2 + Yd21^2) + 
    2*(Yd00*(md200*Yd01 + md210*Yd11 + md220*Yd21) + 
      Yd10*(md201*Yd01 + md211*Yd11 + md221*Yd21) + 
      Yd20*(md202*Yd01 + md212*Yd11 + md222*Yd21)) + 
    Yd20*(mq210*Yd20 + mq211*Yd21 + mq212*Yd22) + 
    mq220*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + 
    Yu00*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
    Yu10*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
    2*mHu2*(Yu00*Yu01 + Yu10*Yu11 + Yu20*Yu21) + 
    mq200*(Yu00*Yu01 + Yu10*Yu11 + Yu20*Yu21) + 
    mq210*(Yu01^2 + Yu11^2 + Yu21^2) + 
    2*(Yu00*(mu200*Yu01 + mu210*Yu11 + mu220*Yu21) + 
      Yu10*(mu201*Yu01 + mu211*Yu11 + mu221*Yu21) + 
      Yu20*(mu202*Yu01 + mu212*Yu11 + mu222*Yu21)) + 
    Yu20*(mq210*Yu20 + mq211*Yu21 + mq212*Yu22) + 
    mq220*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22)) + 
  ((-3*g1^2)/5 + 12*gN^2*QH1p*QQp + 6*(Yd01^2 + Yd11^2 + Yd21^2))*
   ((-2*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 - 6*g2^2*MassWB^2 + 
    (g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QQp^2 + 2*gN^2*QQp*(3*(md200 + md211 + md222)*Qdp + 
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
  6*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22)*
   (2*(TYd01*TYd02 + TYd11*TYd12 + TYd21*TYd22) + 
    2*(TYu01*TYu02 + TYu11*TYu12 + TYu21*TYu22) + 
    Yd02*(mq210*Yd00 + mq211*Yd01 + mq212*Yd02) + 
    Yd12*(mq210*Yd10 + mq211*Yd11 + mq212*Yd12) + 
    mq202*(Yd00*Yd01 + Yd10*Yd11 + Yd20*Yd21) + 
    mq212*(Yd01^2 + Yd11^2 + Yd21^2) + Yd22*(mq210*Yd20 + mq211*Yd21 + 
      mq212*Yd22) + 2*mHd2*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + 
    mq222*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + 
    2*(Yd02*(md200*Yd01 + md210*Yd11 + md220*Yd21) + 
      Yd12*(md201*Yd01 + md211*Yd11 + md221*Yd21) + 
      (md202*Yd01 + md212*Yd11 + md222*Yd21)*Yd22) + 
    Yu02*(mq210*Yu00 + mq211*Yu01 + mq212*Yu02) + 
    Yu12*(mq210*Yu10 + mq211*Yu11 + mq212*Yu12) + 
    mq202*(Yu00*Yu01 + Yu10*Yu11 + Yu20*Yu21) + 
    mq212*(Yu01^2 + Yu11^2 + Yu21^2) + Yu22*(mq210*Yu20 + mq211*Yu21 + 
      mq212*Yu22) + 2*mHu2*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 
    mq222*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 
    2*(Yu02*(mu200*Yu01 + mu210*Yu11 + mu220*Yu21) + 
      Yu12*(mu201*Yu01 + mu211*Yu11 + mu221*Yu21) + 
      (mu202*Yu01 + mu212*Yu11 + mu222*Yu21)*Yu22)) + 
  (4*Lambdax*mHd2 + 4*Lambdax*mHu2 + 4*Lambdax*ms2)*
   ((-3*g1^2*Lambdax)/5 - 3*g2^2*Lambdax + 
    3*(Kappa00^2 + Kappa01^2 + Kappa02^2 + Kappa10^2 + Kappa11^2 + 
      Kappa12^2 + Kappa20^2 + Kappa21^2 + Kappa22^2)*Lambdax + 
    2*(Lambda1200^2 + Lambda1201^2 + Lambda1210^2 + Lambda1211^2)*Lambdax + 
    4*Lambdax^3 - 2*gN^2*Lambdax*QH1p^2 - 2*gN^2*Lambdax*QH2p^2 - 
    2*gN^2*Lambdax*QSp^2 + 3*Lambdax*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    Lambdax*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + 3*Lambdax*(Yu00^2 + Yu01^2 + Yu02^2 + Yu10^2 + 
      Yu11^2 + Yu12^2 + Yu20^2 + Yu21^2 + Yu22^2)) + 
  6*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22)*
   (2*(TYd00*TYd02 + TYd10*TYd12 + TYd20*TYd22) + 
    2*(TYu00*TYu02 + TYu10*TYu12 + TYu20*TYu22) + 
    Yd00*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
    Yd10*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
    Yd20*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22) + 
    2*mHd2*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + 
    mq200*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + 
    mq210*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + 
    mq220*(Yd02^2 + Yd12^2 + Yd22^2) + 
    2*(Yd00*(md200*Yd02 + md210*Yd12 + md220*Yd22) + 
      Yd10*(md201*Yd02 + md211*Yd12 + md221*Yd22) + 
      Yd20*(md202*Yd02 + md212*Yd12 + md222*Yd22)) + 
    Yu00*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
    Yu10*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
    Yu20*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22) + 
    2*mHu2*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 
    mq200*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 
    mq210*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 
    mq220*(Yu02^2 + Yu12^2 + Yu22^2) + 
    2*(Yu00*(mu200*Yu02 + mu210*Yu12 + mu220*Yu22) + 
      Yu10*(mu201*Yu02 + mu211*Yu12 + mu221*Yu22) + 
      Yu20*(mu202*Yu02 + mu212*Yu12 + mu222*Yu22))) + 
  6*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22)*
   (2*(TYd01*TYd02 + TYd11*TYd12 + TYd21*TYd22) + 
    2*(TYu01*TYu02 + TYu11*TYu12 + TYu21*TYu22) + 
    Yd01*(mq220*Yd00 + mq221*Yd01 + mq222*Yd02) + 
    Yd11*(mq220*Yd10 + mq221*Yd11 + mq222*Yd12) + 
    Yd21*(mq220*Yd20 + mq221*Yd21 + mq222*Yd22) + 
    mq201*(Yd00*Yd02 + Yd10*Yd12 + Yd20*Yd22) + 
    2*mHd2*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + 
    mq211*(Yd01*Yd02 + Yd11*Yd12 + Yd21*Yd22) + 
    mq221*(Yd02^2 + Yd12^2 + Yd22^2) + 
    2*(Yd01*(md200*Yd02 + md210*Yd12 + md220*Yd22) + 
      Yd11*(md201*Yd02 + md211*Yd12 + md221*Yd22) + 
      Yd21*(md202*Yd02 + md212*Yd12 + md222*Yd22)) + 
    Yu01*(mq220*Yu00 + mq221*Yu01 + mq222*Yu02) + 
    Yu11*(mq220*Yu10 + mq221*Yu11 + mq222*Yu12) + 
    Yu21*(mq220*Yu20 + mq221*Yu21 + mq222*Yu22) + 
    mq201*(Yu00*Yu02 + Yu10*Yu12 + Yu20*Yu22) + 
    2*mHu2*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 
    mq211*(Yu01*Yu02 + Yu11*Yu12 + Yu21*Yu22) + 
    mq221*(Yu02^2 + Yu12^2 + Yu22^2) + 
    2*(Yu01*(mu200*Yu02 + mu210*Yu12 + mu220*Yu22) + 
      Yu11*(mu201*Yu02 + mu211*Yu12 + mu221*Yu22) + 
      Yu21*(mu202*Yu02 + mu212*Yu12 + mu222*Yu22))) + 
  ((-3*g1^2)/5 + 12*gN^2*QH1p*QQp + 6*(Yd02^2 + Yd12^2 + Yd22^2))*
   ((-2*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 - 6*g2^2*MassWB^2 + 
    (g1^2*(md200 + md211 + md222 - mDx200 - mDx211 - mDx222 + mDxbar200 + 
       mDxbar211 + mDxbar222 + me200 + me211 + me222 - mH1I200 - mH1I211 + 
       mH2I200 + mH2I211 - mHd2 - mHp2 + mHpbar2 + mHu2 - ml200 - ml211 - 
       ml222 + mq200 + mq211 + mq222 - 2*(mu200 + mu211 + mu222)))/5 - 
    8*gN^2*MassBp^2*QQp^2 + 2*gN^2*QQp*(3*(md200 + md211 + md222)*Qdp + 
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
  ((-3*g1^2)/5 + 2*Lambdax^2 + 4*gN^2*QH1p*QH2p)*
   ((-6*g1^2*MassB^2)/5 - 6*g2^2*MassWB^2 + 2*Lambdax^2*mHd2 + 
    2*Lambdax^2*mHu2 + 2*Lambdax^2*ms2 + 
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
  12*TYd00*((-7*g1^2*TYd00)/15 - 3*g2^2*TYd00 - (16*g3^2*TYd00)/3 + 
    Lambdax^2*TYd00 - 2*gN^2*Qdp^2*TYd00 - 2*gN^2*QH1p^2*TYd00 - 
    2*gN^2*QQp^2*TYd00 + 3*TYd00*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    5*(Yd00*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02) + 
      Yd10*(TYd00*Yd10 + TYd01*Yd11 + TYd02*Yd12) + 
      Yd20*(TYd00*Yd20 + TYd01*Yd21 + TYd02*Yd22)) + 
    4*(TYd00*(Yd00^2 + Yd01^2 + Yd02^2) + TYd10*(Yd00*Yd10 + Yd01*Yd11 + 
        Yd02*Yd12) + TYd20*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    TYd00*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + Yd00*((14*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      6*g2^2*MassWB + 4*gN^2*MassBp*Qdp^2 + 4*gN^2*MassBp*QH1p^2 + 
      4*gN^2*MassBp*QQp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Yu00*(TYd00*Yu00 + TYd01*Yu01 + TYd02*Yu02) + 
    Yu10*(TYd00*Yu10 + TYd01*Yu11 + TYd02*Yu12) + 
    Yu20*(TYd00*Yu20 + TYd01*Yu21 + TYd02*Yu22) + 
    2*(TYu00*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      TYu10*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      TYu20*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22))) + 
  12*TYd01*((-7*g1^2*TYd01)/15 - 3*g2^2*TYd01 - (16*g3^2*TYd01)/3 + 
    Lambdax^2*TYd01 - 2*gN^2*Qdp^2*TYd01 - 2*gN^2*QH1p^2*TYd01 - 
    2*gN^2*QQp^2*TYd01 + 3*TYd01*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    5*(Yd01*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02) + 
      Yd11*(TYd00*Yd10 + TYd01*Yd11 + TYd02*Yd12) + 
      Yd21*(TYd00*Yd20 + TYd01*Yd21 + TYd02*Yd22)) + 
    4*(TYd01*(Yd00^2 + Yd01^2 + Yd02^2) + TYd11*(Yd00*Yd10 + Yd01*Yd11 + 
        Yd02*Yd12) + TYd21*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    TYd01*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + Yd01*((14*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      6*g2^2*MassWB + 4*gN^2*MassBp*Qdp^2 + 4*gN^2*MassBp*QH1p^2 + 
      4*gN^2*MassBp*QQp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Yu01*(TYd00*Yu00 + TYd01*Yu01 + TYd02*Yu02) + 
    Yu11*(TYd00*Yu10 + TYd01*Yu11 + TYd02*Yu12) + 
    Yu21*(TYd00*Yu20 + TYd01*Yu21 + TYd02*Yu22) + 
    2*(TYu01*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      TYu11*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      TYu21*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22))) + 
  12*TYd02*((-7*g1^2*TYd02)/15 - 3*g2^2*TYd02 - (16*g3^2*TYd02)/3 + 
    Lambdax^2*TYd02 - 2*gN^2*Qdp^2*TYd02 - 2*gN^2*QH1p^2*TYd02 - 
    2*gN^2*QQp^2*TYd02 + 3*TYd02*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    5*(Yd02*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02) + 
      Yd12*(TYd00*Yd10 + TYd01*Yd11 + TYd02*Yd12) + 
      Yd22*(TYd00*Yd20 + TYd01*Yd21 + TYd02*Yd22)) + 
    4*(TYd02*(Yd00^2 + Yd01^2 + Yd02^2) + TYd12*(Yd00*Yd10 + Yd01*Yd11 + 
        Yd02*Yd12) + TYd22*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22)) + 
    TYd02*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + Yd02*((14*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      6*g2^2*MassWB + 4*gN^2*MassBp*Qdp^2 + 4*gN^2*MassBp*QH1p^2 + 
      4*gN^2*MassBp*QQp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Yu02*(TYd00*Yu00 + TYd01*Yu01 + TYd02*Yu02) + 
    Yu12*(TYd00*Yu10 + TYd01*Yu11 + TYd02*Yu12) + 
    Yu22*(TYd00*Yu20 + TYd01*Yu21 + TYd02*Yu22) + 
    2*(TYu02*(Yd00*Yu00 + Yd01*Yu01 + Yd02*Yu02) + 
      TYu12*(Yd00*Yu10 + Yd01*Yu11 + Yd02*Yu12) + 
      TYu22*(Yd00*Yu20 + Yd01*Yu21 + Yd02*Yu22))) + 
  12*TYd10*((-7*g1^2*TYd10)/15 - 3*g2^2*TYd10 - (16*g3^2*TYd10)/3 + 
    Lambdax^2*TYd10 - 2*gN^2*Qdp^2*TYd10 - 2*gN^2*QH1p^2*TYd10 - 
    2*gN^2*QQp^2*TYd10 + 3*TYd10*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    5*(Yd00*(TYd10*Yd00 + TYd11*Yd01 + TYd12*Yd02) + 
      Yd10*(TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12) + 
      Yd20*(TYd10*Yd20 + TYd11*Yd21 + TYd12*Yd22)) + 
    4*(TYd00*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      TYd10*(Yd10^2 + Yd11^2 + Yd12^2) + TYd20*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + TYd10*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yd10*((14*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 6*g2^2*MassWB + 
      4*gN^2*MassBp*Qdp^2 + 4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QQp^2 + 
      2*Lambdax*TLambdax + 6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + 
        TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + 
        TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
        TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Yu00*(TYd10*Yu00 + TYd11*Yu01 + TYd12*Yu02) + 
    Yu10*(TYd10*Yu10 + TYd11*Yu11 + TYd12*Yu12) + 
    Yu20*(TYd10*Yu20 + TYd11*Yu21 + TYd12*Yu22) + 
    2*(TYu00*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      TYu10*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      TYu20*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22))) + 
  12*TYd11*((-7*g1^2*TYd11)/15 - 3*g2^2*TYd11 - (16*g3^2*TYd11)/3 + 
    Lambdax^2*TYd11 - 2*gN^2*Qdp^2*TYd11 - 2*gN^2*QH1p^2*TYd11 - 
    2*gN^2*QQp^2*TYd11 + 3*TYd11*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    5*(Yd01*(TYd10*Yd00 + TYd11*Yd01 + TYd12*Yd02) + 
      Yd11*(TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12) + 
      Yd21*(TYd10*Yd20 + TYd11*Yd21 + TYd12*Yd22)) + 
    4*(TYd01*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      TYd11*(Yd10^2 + Yd11^2 + Yd12^2) + TYd21*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + TYd11*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yd11*((14*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 6*g2^2*MassWB + 
      4*gN^2*MassBp*Qdp^2 + 4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QQp^2 + 
      2*Lambdax*TLambdax + 6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + 
        TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + 
        TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
        TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Yu01*(TYd10*Yu00 + TYd11*Yu01 + TYd12*Yu02) + 
    Yu11*(TYd10*Yu10 + TYd11*Yu11 + TYd12*Yu12) + 
    Yu21*(TYd10*Yu20 + TYd11*Yu21 + TYd12*Yu22) + 
    2*(TYu01*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      TYu11*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      TYu21*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22))) + 
  12*TYd12*((-7*g1^2*TYd12)/15 - 3*g2^2*TYd12 - (16*g3^2*TYd12)/3 + 
    Lambdax^2*TYd12 - 2*gN^2*Qdp^2*TYd12 - 2*gN^2*QH1p^2*TYd12 - 
    2*gN^2*QQp^2*TYd12 + 3*TYd12*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    5*(Yd02*(TYd10*Yd00 + TYd11*Yd01 + TYd12*Yd02) + 
      Yd12*(TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12) + 
      Yd22*(TYd10*Yd20 + TYd11*Yd21 + TYd12*Yd22)) + 
    4*(TYd02*(Yd00*Yd10 + Yd01*Yd11 + Yd02*Yd12) + 
      TYd12*(Yd10^2 + Yd11^2 + Yd12^2) + TYd22*(Yd10*Yd20 + Yd11*Yd21 + 
        Yd12*Yd22)) + TYd12*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + 
      Ye12^2 + Ye20^2 + Ye21^2 + Ye22^2) + 
    Yd12*((14*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 6*g2^2*MassWB + 
      4*gN^2*MassBp*Qdp^2 + 4*gN^2*MassBp*QH1p^2 + 4*gN^2*MassBp*QQp^2 + 
      2*Lambdax*TLambdax + 6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + 
        TYd10*Yd10 + TYd11*Yd11 + TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + 
        TYd22*Yd22) + 2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + 
        TYe11*Ye11 + TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Yu02*(TYd10*Yu00 + TYd11*Yu01 + TYd12*Yu02) + 
    Yu12*(TYd10*Yu10 + TYd11*Yu11 + TYd12*Yu12) + 
    Yu22*(TYd10*Yu20 + TYd11*Yu21 + TYd12*Yu22) + 
    2*(TYu02*(Yd10*Yu00 + Yd11*Yu01 + Yd12*Yu02) + 
      TYu12*(Yd10*Yu10 + Yd11*Yu11 + Yd12*Yu12) + 
      TYu22*(Yd10*Yu20 + Yd11*Yu21 + Yd12*Yu22))) + 
  12*TYd20*((-7*g1^2*TYd20)/15 - 3*g2^2*TYd20 - (16*g3^2*TYd20)/3 + 
    Lambdax^2*TYd20 - 2*gN^2*Qdp^2*TYd20 - 2*gN^2*QH1p^2*TYd20 - 
    2*gN^2*QQp^2*TYd20 + 3*TYd20*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    5*(Yd00*(TYd20*Yd00 + TYd21*Yd01 + TYd22*Yd02) + 
      Yd10*(TYd20*Yd10 + TYd21*Yd11 + TYd22*Yd12) + 
      Yd20*(TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)) + 
    4*(TYd00*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      TYd10*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      TYd20*(Yd20^2 + Yd21^2 + Yd22^2)) + 
    TYd20*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + Yd20*((14*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      6*g2^2*MassWB + 4*gN^2*MassBp*Qdp^2 + 4*gN^2*MassBp*QH1p^2 + 
      4*gN^2*MassBp*QQp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Yu00*(TYd20*Yu00 + TYd21*Yu01 + TYd22*Yu02) + 
    Yu10*(TYd20*Yu10 + TYd21*Yu11 + TYd22*Yu12) + 
    Yu20*(TYd20*Yu20 + TYd21*Yu21 + TYd22*Yu22) + 
    2*(TYu00*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      TYu10*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      TYu20*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) + 
  12*TYd21*((-7*g1^2*TYd21)/15 - 3*g2^2*TYd21 - (16*g3^2*TYd21)/3 + 
    Lambdax^2*TYd21 - 2*gN^2*Qdp^2*TYd21 - 2*gN^2*QH1p^2*TYd21 - 
    2*gN^2*QQp^2*TYd21 + 3*TYd21*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    5*(Yd01*(TYd20*Yd00 + TYd21*Yd01 + TYd22*Yd02) + 
      Yd11*(TYd20*Yd10 + TYd21*Yd11 + TYd22*Yd12) + 
      Yd21*(TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)) + 
    4*(TYd01*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      TYd11*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      TYd21*(Yd20^2 + Yd21^2 + Yd22^2)) + 
    TYd21*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + Yd21*((14*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      6*g2^2*MassWB + 4*gN^2*MassBp*Qdp^2 + 4*gN^2*MassBp*QH1p^2 + 
      4*gN^2*MassBp*QQp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Yu01*(TYd20*Yu00 + TYd21*Yu01 + TYd22*Yu02) + 
    Yu11*(TYd20*Yu10 + TYd21*Yu11 + TYd22*Yu12) + 
    Yu21*(TYd20*Yu20 + TYd21*Yu21 + TYd22*Yu22) + 
    2*(TYu01*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      TYu11*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      TYu21*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) + 
  12*TYd22*((-7*g1^2*TYd22)/15 - 3*g2^2*TYd22 - (16*g3^2*TYd22)/3 + 
    Lambdax^2*TYd22 - 2*gN^2*Qdp^2*TYd22 - 2*gN^2*QH1p^2*TYd22 - 
    2*gN^2*QQp^2*TYd22 + 3*TYd22*(Yd00^2 + Yd01^2 + Yd02^2 + Yd10^2 + 
      Yd11^2 + Yd12^2 + Yd20^2 + Yd21^2 + Yd22^2) + 
    5*(Yd02*(TYd20*Yd00 + TYd21*Yd01 + TYd22*Yd02) + 
      Yd12*(TYd20*Yd10 + TYd21*Yd11 + TYd22*Yd12) + 
      Yd22*(TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22)) + 
    4*(TYd02*(Yd00*Yd20 + Yd01*Yd21 + Yd02*Yd22) + 
      TYd12*(Yd10*Yd20 + Yd11*Yd21 + Yd12*Yd22) + 
      TYd22*(Yd20^2 + Yd21^2 + Yd22^2)) + 
    TYd22*(Ye00^2 + Ye01^2 + Ye02^2 + Ye10^2 + Ye11^2 + Ye12^2 + Ye20^2 + 
      Ye21^2 + Ye22^2) + Yd22*((14*g1^2*MassB)/15 + (32*g3^2*MassG)/3 + 
      6*g2^2*MassWB + 4*gN^2*MassBp*Qdp^2 + 4*gN^2*MassBp*QH1p^2 + 
      4*gN^2*MassBp*QQp^2 + 2*Lambdax*TLambdax + 
      6*(TYd00*Yd00 + TYd01*Yd01 + TYd02*Yd02 + TYd10*Yd10 + TYd11*Yd11 + 
        TYd12*Yd12 + TYd20*Yd20 + TYd21*Yd21 + TYd22*Yd22) + 
      2*(TYe00*Ye00 + TYe01*Ye01 + TYe02*Ye02 + TYe10*Ye10 + TYe11*Ye11 + 
        TYe12*Ye12 + TYe20*Ye20 + TYe21*Ye21 + TYe22*Ye22)) + 
    Yu02*(TYd20*Yu00 + TYd21*Yu01 + TYd22*Yu02) + 
    Yu12*(TYd20*Yu10 + TYd21*Yu11 + TYd22*Yu12) + 
    Yu22*(TYd20*Yu20 + TYd21*Yu21 + TYd22*Yu22) + 
    2*(TYu02*(Yd20*Yu00 + Yd21*Yu01 + Yd22*Yu02) + 
      TYu12*(Yd20*Yu10 + Yd21*Yu11 + Yd22*Yu12) + 
      TYu22*(Yd20*Yu20 + Yd21*Yu21 + Yd22*Yu22))) + 
  ((6*g1^2)/5 + 6*gN^2*QH1p*Qup)*((-32*g1^2*MassB^2)/15 - 
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
    4*(TYu00^2 + TYu01^2 + TYu02^2) + 4*mHu2*(Yu00^2 + Yu01^2 + Yu02^2) + 
    4*(Yu00*(mq200*Yu00 + mq210*Yu01 + mq220*Yu02) + 
      Yu01*(mq201*Yu00 + mq211*Yu01 + mq221*Yu02) + 
      Yu02*(mq202*Yu00 + mq212*Yu01 + mq222*Yu02)) + 
    2*(Yu00*(mu200*Yu00 + mu201*Yu10 + mu202*Yu20) + 
      Yu01*(mu200*Yu01 + mu201*Yu11 + mu202*Yu21) + 
      Yu02*(mu200*Yu02 + mu201*Yu12 + mu202*Yu22)) + 
    2*(mu200*(Yu00^2 + Yu01^2 + Yu02^2) + mu210*(Yu00*Yu10 + Yu01*Yu11 + 
        Yu02*Yu12) + mu220*(Yu00*Yu20 + Yu01*Yu21 + Yu02*Yu22))) + 
  ((6*g1^2)/5 + 6*gN^2*QH1p*Qup)*((-32*g1^2*MassB^2)/15 - 
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
    4*(TYu10^2 + TYu11^2 + TYu12^2) + 4*mHu2*(Yu10^2 + Yu11^2 + Yu12^2) + 
    4*(Yu10*(mq200*Yu10 + mq210*Yu11 + mq220*Yu12) + 
      Yu11*(mq201*Yu10 + mq211*Yu11 + mq221*Yu12) + 
      Yu12*(mq202*Yu10 + mq212*Yu11 + mq222*Yu12)) + 
    2*(Yu10*(mu210*Yu00 + mu211*Yu10 + mu212*Yu20) + 
      Yu11*(mu210*Yu01 + mu211*Yu11 + mu212*Yu21) + 
      Yu12*(mu210*Yu02 + mu211*Yu12 + mu212*Yu22)) + 
    2*(mu201*(Yu00*Yu10 + Yu01*Yu11 + Yu02*Yu12) + 
      mu211*(Yu10^2 + Yu11^2 + Yu12^2) + mu221*(Yu10*Yu20 + Yu11*Yu21 + 
        Yu12*Yu22))) + ((6*g1^2)/5 + 6*gN^2*QH1p*Qup)*
   ((-32*g1^2*MassB^2)/15 - (32*g3^2*MassG^2)/3 - 
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
  4*TLambdax*((6*g1^2*Lambdax*MassB)/5 + 6*g2^2*Lambdax*MassWB + 
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
