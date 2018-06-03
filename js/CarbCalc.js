// initially, code to re-create pH isopleths on WQ map

const CarbCalc = {

// # nb: pH(NBS) = pH(Free) -
// #               log10(ahFreeToSwsFactor(S, T(Kelvin), 0)) -
// #               log10(ahSwsToNbsFactor())
    phNbsToPhFree: (ph, sal, temp, p) => {
      return(ph +
               Math.log10(ahSwsToNbsFactor(sal, temp, p)) +
               Math.log10(ahFreeToSwsFactor(sal, temp, p))
             );
    },

    ahSwsToNbsFactor: (sal, temp, p) => {
      return(calcProtonActivityCoeffZg(temp, sal, 0.0) / ahMolalToMolinforSalinity(sal));
    },

    calcProtonActivityCoeffZg: (temp, sal, p) => {

      let rootGamma = Math.sqrt(calcIonicStrength(sal));

      let myHActivityCoeff = 1820000.0 * Math.pow((79 * temp), -1.5);

      myHActivityCoeff = myHActivityCoeff * ((rootGamma / (1 + rootGamma)) - 0.2 * calcIonicStrength(sal));

      myHActivityCoeff = Math.pow(10, -myHActivityCoeff);

      return(myHActivityCoeff);
    },

    calcIonicStrength: (sal) => {
      let myIS = 19.924 * sal / (1000.0 - 1.005 * sal); // mole/kg-H2O (molal)
      return(myIS);
    },

    ahMolalToMolinforSalinity: (sal) => {
      return(1.0 - 0.001005 * sal);
    },


    ahFreeToSwsFactor: (sal, temp, p) => {
      return(1 + (calcTS(sal) / calcKsDickson(temp, sal, 0.0)) +
               (calcTF(sal) / calcKfDickson(temp, sal, 0.0)));
    },

    calcTF: (sal) => {
      return(0.0000019522 * sal);
    },

    calcTS: (sal) => {
      return(0.0008067267 * sal);
    },

    // # pH scale: FREE
    // # concentration scale: mole/kg-H2O ?? -> molin...?
    calcKfDickson: (temp, sal, p) => {

      let fluorFactor1 = 1590.2 / temp;
      let fluorFactor2 = -12.641;
      let fluorFactor3 = 1.525 * Math.sqrt(calcIonicStrength(sal));

      // # ** nb: molalToMolin factor in sqrt()...
      // #    double fluorFactor3 = 1.525 * sqrt([self calcIonicStrength:s] * (1.0 - 0.001005 * s));

      let molal2molin = Math.log(1.0 - 0.001005 * sal);

      let KF = fluorFactor1 + fluorFactor2 + fluorFactor3 + molal2molin;

      return(Math.exp(KF))
    },

    // # bisulfate dissociation
    // # Dickson (1990); DOE (1994), ch. 5 p. 13; Z & W-G (2001) p. 260
    // # pH scale: Free
    // # concentration scale: mol/kg-H2O, CONVERTED TO AND RETURNED AS mol/kg-soln (molin)
    // # **** nb **** if called from within pH conversion, T already C -> K *******
    // # **** nb **** else if called from getKS(), must add 273.15 in call ****************
    calcKsDickson: (temp,sal,p) => {

      // # ** NB: ionic strength calc now returns molal, so
      // # **     change to molin (mol/kg-soln) here for this calc
      // #    double myIS = [self calcIonicStrength:s] * (1.0 - 0.001005 * s);
      let myIS = calcIonicStrength(sal);

      let sulfFactor1 = 141.328 - 4276.1 / temp;
      let sulfFactor2 = -23.093 * Math.log(temp);
      let sulfFactor3 = ( 324.57 - 47.986 * Math.log(temp) - 13856 / temp) * Math.sqrt(myIS)
      let sulfFactor4 = (-771.54 + 114.723 * Math.log(temp) + 35474 / temp) * myIS
      let sulfFactor5 = -2698.0 * Math.pow(myIS, 1.5) / temp
      let sulfFactor6 =  1776.0 * myIS * myIS / temp

      molal2molin = Math.log(1.0 - 0.001005 * sal)

     // # ** in mol/kg-soln (MOLIN)
      let KS = sulfFactor1 + sulfFactor2 + sulfFactor3 + sulfFactor4 + sulfFactor5 + sulfFactor6 + molal2molin;

     // # ** in mol/kg-H2O (MOLAL)
     // #	double KS = sulfFactor1 + sulfFactor2 + sulfFactor3 + sulfFactor4 + sulfFactor5 + sulfFactor6;

      return(Math.exp(KS))
    },


    phLineIntercept: (temp, sal, ph) => {
      // console.log(`calcHydroxide = ${calcHydroxide(ph, temp, sal)}`);
      // console.log(`calcHydronium = ${calcHydronium(ph)}`);
      // console.log(`   calcBorate = ${calcBorate(ph, temp, sal)}`);
      return(calcHydroxide(ph, temp, sal) - calcHydronium(ph) + calcBorate(ph, temp, sal));
    },

    calcHydroxide: (ph, t, sal) => {
      // if(sal === null)
      //   return();
      let kWToTheTen = Math.log10(calcKWMehrbach(t, sal));

      return(Math.pow(10, (kWToTheTen + ph)));
    },

    calcKWMehrbach: (temp, sal) => {
      // if(sal === null)
      //   return()

      let expSum = 148.9652;     // # nb: "148.965 02" in Zeebe & Wolf-Gladrow code
      expSum = expSum - 13847.26 / temp;
      expSum = expSum - 23.6521 * Math.log(temp);
      expSum = expSum + (-5.977 + (118.67 / temp) + 1.0495 * Math.log(temp)) * Math.sqrt(sal);
      expSum = expSum - 0.01615 * sal;

      let KW = Math.exp(expSum);   // # still on TOTAL scale
      let pKW = -Math.log10(KW);   // # still on TOTAL scale

      // # ** nb: convert to FREE pH scale, as per AquaEnv
      pKW = pKW + Math.log10(ahFreeToTotFactor(sal, temp, 0.0));

      return(Math.pow(10, -pKW));
    },

    ahFreeToTotFactor: (sal, temp, p) => {
      return(1 + (calcTS(sal) / calcKsDickson(temp, sal, 0.0)));
    },

    calcHydronium: (ph) => {
      return(Math.pow(10, -ph));
    },

    phLineSlope: (temp, sal, ph) => {
      return(alphaOne(temp, sal, ph) + 2 * alphaTwo(temp, sal, ph));
    },

    // # Borate ----
    // # ** FREE scale?
    calcBorate: (ph, t, sal) => {

      let concB = calcBorateConcOfSalinity(sal);

      let myKB = calcBorateFactor(t, sal);

      let borate = myKB * concB / (myKB + Math.pow(10, -ph));

      return(borate);
    },


    calcBorateConcOfSalinity: (sal) => {

      let concB = 0.000232 * sal / (10.811 * 1.80655);

      return(concB);
    },

    calcBorateFactor: (t, sal) => {
      // if(sal === null)
      //   return();

      let A = 148.0248 + 137.1942 * Math.sqrt(sal) + 1.62142 * sal;
      let B = -8966.90 - 2890.53 * Math.sqrt(sal) - 77.942 * sal + 1.728 * Math.pow(sal, 1.5) - 0.0996 * sal * sal;
      let C = -24.4344 - 25.085 * Math.sqrt(sal) - 0.2474 * sal;
      let D = 0.053105 * Math.sqrt(sal);

      let K_BOH3 = Math.exp(A + B/t + C * Math.log(t) + D * t);

      let ans = Math.pow(10, -(-Math.log10(K_BOH3) + Math.log10(ahFreeToTotFactor(sal, t, 0))));

      return(ans);
    },


    alphaOne: (temp, sal, ph) => {

      // # nb: define p LOCALLY until incorporate in calcs
      let p = 0.0;

      let h = calcHydronium(ph);

      let k1 = getK1(temp, sal, 0);

      let k2 = getK2(temp, sal, 0);

      let numerator = h * k1;

      return (numerator / calcAlphaDenom(h, k1, k2));
    },

    alphaTwo: function(temp, sal, ph) {

      // # nb: define p LOCALLY until incorporate in calcs
      let p = 0.0;

      let h = calcHydronium(ph);

      let k1 = getK1(temp, sal, 0);

      let k2 = getK2(temp, sal, 0);

      let numerator = k1 * k2;

      return (numerator / calcAlphaDenom(h, k1, k2));
    },

    calcAlphaDenom: (h, k1, k2) => {
      return(h * h + k1 * h + k1 * k2);
    },

    // # T- and S-dependent K1 from Millero et al. (2006)
    // # ** pH scale: SWS for calculation, FREE returned
    getK1: (temp, sal, p) => {
      // if(is.null(sal))
      //   return()

      sqrtS = Math.sqrt(sal);
      lnT   = Math.log(temp);

      let pK1z = -126.34048 + (6320.813 / temp) + 19.568224 * lnT;  // # was "19.56822.."
      let A =   13.4191 * sqrtS + 0.0331 * sal - 0.0000533 * sal * sal;
      let B = -530.123 * sqrtS - 6.103 * sal;
      let C =   -2.06950 * sqrtS;
      let pK1 = pK1z + A + (B / temp) + C * lnT;

      // # ** nb: THIS conversion puts it on the FREE pH scale *from* the SWS scale, as per AquaEnv
      // # **     factor for sws2free = 1.0 / ahFreeToSwsFactorForSalinity:temp:pressure:
      // # **     when dealing with -log10, -log10(1 / ah) = -(-log10(ah)) = +log10(ahFreeToSWSFactor...)
      pK1 = pK1 + Math.log10(ahFreeToSwsFactor(sal, temp, 0.0));

      // # ?? nb: no concentration scale conversion needed, as both on molinity [sic]

      return(Math.pow(10, -pK1));  // # K1
    },


    // # T- and S-dependent K2 from Millero et al. (2006)
    // # ** pH scale: SWS for calculation
    // # ** Return: K2 -- not pK2 -- on FREE pH scale

    // # ----> nb: Millero (2010) slightly changes some coefficients <----

    getK2: (temp, sal, p) => {

      // if(is.null(sal))
      //   return()

      let sqrtS = Math.sqrt(sal);
      let lnT   = Math.log(temp);

      let pK2z = -90.18333 + (5143.692 / temp) + 14.613358 * lnT;

      let A = 21.0894 * sqrtS + 0.1248 * sal - 0.0003687 * sal * sal; // # was "21.08945"
      let B = -772.483 * sqrtS - 20.051 * sal;
      let C = -3.3336 * sqrtS; //  # was "-3.3336" or "3.32254"

      // # ** nb: THIS pK2 is on the SWS pH scale
      let pK2 = pK2z + A + (B / temp) + C * lnT;

      pK2 = pK2 + Math.log10(ahFreeToSwsFactor(sal, temp, 0.0));

      return(Math.pow(10, -pK2));
    },


    calcDicOfAlk: (alk, ph, temp, sal) => {

      let m = phLineSlope(temp, sal, ph);

      let dic = (alk - calcHydroxide(ph, temp, sal) -
                      calcBorate(ph, temp, sal) +
                      calcHydronium(ph)) / m;
      return(dic);
    }

} // END CarbCalc

// export default CarbCalc;
