  // Isotopes

  std::map<std::string,int> Zs;
  std::map<std::string,int> Ns;
  std::map<std::string,double> Deltas;
  std::map<std::string,int> Zdiffs;
  std::map<std::string,int> Ndiffs;

  Zs["H1"] = 1;
  Zs["H2"] = 1;

  Zs["He3"] = 2;
  Zs["He4"] = 2;

  Zs["C12"] = 6;
  Zs["C13"] = 6;

  Zs["F19"] = 9;

  Zs["Ne20"]= 10;
  Zs["Ne21"]= 10;
  Zs["Ne22"]= 10;

  Zs["Na23"]= 11;

  Zs["Si28"]= 14;
  Zs["Si29"]= 14;
  Zs["Si30"]= 14;

  Zs["S32"] = 16;
  Zs["S33"] = 16;
  Zs["S34"] = 16;
  Zs["S36"] = 16;

  Zs["Ar36"]= 18;
  Zs["Ar38"]= 18;
  Zs["Ar40"]= 18;

  Zs["Ge70"]= 32;
  Zs["Ge72"]= 32;
  Zs["Ge73"]= 32;
  Zs["Ge74"]= 32;
  Zs["Ge76"]= 32;

  Zs["I127"]= 53;

  Zs["Xe124"]= 54;
  Zs["Xe126"]= 54;
  Zs["Xe128"]= 54;
  Zs["Xe129"]= 54;
  Zs["Xe130"]= 54;
  Zs["Xe131"]= 54;
  Zs["Xe132"]= 54;
  Zs["Xe134"]= 54;
  Zs["Xe136"]= 54;

  Zs["Cs133"]= 55;

  Ns["H1"] = 0;
  Ns["H2"] = 1;

  Ns["He3"] = 1;
  Ns["He4"] = 2;

  Ns["C12"] = 6;
  Ns["C13"] = 7;

  Ns["F19"] = 10;

  Ns["Ne20"]= 10;
  Ns["Ne21"]= 11;
  Ns["Ne22"]= 12;

  Ns["Na23"]= 12;

  Ns["Si28"]= 14;
  Ns["Si29"]= 15;
  Ns["Si30"]= 16;

  Ns["S32"] = 16;
  Ns["S33"] = 17;
  Ns["S34"] = 18;
  Ns["S36"] = 20;

  Ns["Ar36"]= 18;
  Ns["Ar38"]= 20;
  Ns["Ar40"]= 22;

  Ns["Ge70"]= 38;
  Ns["Ge72"]= 40;
  Ns["Ge73"]= 41;
  Ns["Ge74"]= 42;
  Ns["Ge76"]= 44;

  Ns["I127"]= 74;

  Ns["Xe124"]= 70;
  Ns["Xe126"]= 72;
  Ns["Xe128"]= 74;
  Ns["Xe129"]= 75;
  Ns["Xe130"]= 76;
  Ns["Xe131"]= 77;
  Ns["Xe132"]= 78;
  Ns["Xe134"]= 80;
  Ns["Xe136"]= 82;

  Ns["Cs133"]= 78;


  Deltas["H1"] = 7.2889;
  Deltas["H2"] = 13.1357;

  Deltas["He3"] = 14.9312;
  Deltas["He4"] = 2.4249;

  Deltas["C12"] = 0;
  Deltas["C13"] = 3.125;

  Deltas["F19"] = -1.4874;

  Deltas["Ne20"]= -7.02;
  Deltas["Ne21"]= -5.731;
  Deltas["Ne22"]= -8.024;

  Deltas["Na23"]= -9.530;


  Deltas["Si28"]= -21.43;
  Deltas["Si29"]= -21.895;
  Deltas["Si30"]= -24.432;

  Deltas["S32"] = -26.0155;
  Deltas["S33"] = -26.5858;
  Deltas["S34"] = -29.9316;
  Deltas["S36"] = -30.6641;


  Deltas["Ar36"]= -30.231;
  Deltas["Ar38"]= -34.714;
  Deltas["Ar40"]= -35.040;

  Deltas["Ge70"]= -70.561;
  Deltas["Ge72"]= -72.585;
  Deltas["Ge73"]= -71.297;
  Deltas["Ge74"]= -73.422;
  Deltas["Ge76"]= -73.212;

  Deltas["I127"]= -88.984;

  Deltas["Xe124"]= -87.661;
  Deltas["Xe126"]= -89.146;
  Deltas["Xe128"]= -89.860;
  Deltas["Xe129"]= -88.696;
  Deltas["Xe130"]= -89.880;
  Deltas["Xe131"]= -88.413;
  Deltas["Xe132"]= -89.279;
  Deltas["Xe134"]= -88.124;
  Deltas["Xe136"]= -86.429;

  Deltas["Cs133"]= -88.070;


  Zdiffs["H1"] = 1;
  Zdiffs["H2"] = 1;

  Zdiffs["He3"] = 0;
  Zdiffs["He4"] = 0;

  Zdiffs["C12"] = 0;
  Zdiffs["C13"] = 0;

  Zdiffs["F19"] = 1;

  Zdiffs["Ne20"]= 0;
  Zdiffs["Ne21"]= 0;
  Zdiffs["Ne22"]= 0;

  Zdiffs["Na23"]= 1;

  Zdiffs["Si28"]= 0;
  Zdiffs["Si29"]= 0;
  Zdiffs["Si30"]= 0;

  Zdiffs["S32"] = 0;
  Zdiffs["S33"] = 0;
  Zdiffs["S34"] = 0;
  Zdiffs["S36"] = 0;

  Zdiffs["Ar36"]= 0;
  Zdiffs["Ar38"]= 0;
  Zdiffs["Ar40"]= 0;

  Zdiffs["Ge70"]= 0;
  Zdiffs["Ge72"]= 0;
  Zdiffs["Ge73"]= 0;
  Zdiffs["Ge74"]= 0;
  Zdiffs["Ge76"]= 0;

  Zdiffs["I127"]= 1;

  Zdiffs["Xe124"]= 0;
  Zdiffs["Xe126"]= 0;
  Zdiffs["Xe128"]= 0;
  Zdiffs["Xe129"]= 0;
  Zdiffs["Xe130"]= 0;
  Zdiffs["Xe131"]= 0;
  Zdiffs["Xe132"]= 0;
  Zdiffs["Xe134"]= 0;
  Zdiffs["Xe136"]= 0;

  Zdiffs["Cs133"]= 1;

  Ndiffs["H1"] = 0;
  Ndiffs["H2"] = 1;

  Ndiffs["He3"] = 1;
  Ndiffs["He4"] = 0;

  Ndiffs["C12"] = 0;
  Ndiffs["C13"] = 1;

  Zdiffs["F19"] = 0;

  Ndiffs["Ne20"]= 0;
  Ndiffs["Ne21"]= 1;
  Ndiffs["Ne22"]= 0;

  Ndiffs["Na23"]= 0;

  Ndiffs["Si28"]= 0;
  Ndiffs["Si29"]= 1;
  Ndiffs["Si30"]= 0;

  Ndiffs["S32"] = 0;
  Ndiffs["S33"] = 1;
  Ndiffs["S34"] = 0;
  Ndiffs["S36"] = 0;


  Ndiffs["Ar36"]= 0;
  Ndiffs["Ar38"]= 0;
  Ndiffs["Ar40"]= 0;

  Ndiffs["Ge70"]= 0;
  Ndiffs["Ge72"]= 0;
  Ndiffs["Ge73"]= 1;
  Ndiffs["Ge74"]= 0;
  Ndiffs["Ge76"]= 0;

  Ndiffs["I127"]= 0;

  Ndiffs["Xe124"]= 0;
  Ndiffs["Xe126"]= 0;
  Ndiffs["Xe128"]= 0;
  Ndiffs["Xe129"]= 1;
  Ndiffs["Xe130"]= 0;
  Ndiffs["Xe131"]= 1;
  Ndiffs["Xe132"]= 0;
  Ndiffs["Xe134"]= 0;
  Ndiffs["Xe136"]= 0;

  Ndiffs["Cs133"]= 0;
