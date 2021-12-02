#include <stdio.h>
#include <math.h>
int main()
{
    float Lva, vat2,va,t2, p, dp, Omg4,eta3,kc;
    float epsb2,epse1;
    float m,mv2,v9,Vdn;
    float G0,G1,G2,G3;
    float gm1,gm2,gm3,gm4,gm5,gm6,gm7,gm8,gm9,gm10;
    int i;
    float beta,dtm;
    float p9,p18,Gms1,Gms2,Gms3,Gms4,Gms5,Gms6,Gms7,Gms8,Gms9,Gms10,Gms11,GP4,GP5;

    printf("Radio Emission from Outflow-Cloud Interaction and Its Constraint on TDE outflow.\n");
    printf("Ref: Mou et al. 2021, arXiv: 2108.11296 \n");
    printf("Purpose: Calculating the physics of TDE outflow by setting up \n");
    printf("         cloud parameters (assumed) and \n");
    printf("         radio parameters (from observations). \n \n");
    printf("Suggestion: you can test it against Table 2 in the above paper before use. \n");
    printf("Caution: by default, epsilon_{B,-2}=Bar{epsilon}_{e,-1}= 1 . \n");
    printf("If you want to test other values of epsilon_{B,-2} and Bar{epsilon}_{e,-1}, \n"); 
    printf("please modify them in this code. epsilon_{B,-2} is epsb2, Bar{epsilon}_{e,-1} is epse1 \n \n");
    epsb2=1.0;
    epse1=1.0;
/* Fiducial values as references */
    Omg4=4.0/4.0; /* Omega_4 */
    eta3=1.0;
    kc=1.0 ;
    
    printf("1. Setting up cloud parameters \n");
    printf("Enter Omega_4 (Omega/4sr): \n");
    scanf("%f", &Omg4);
    if(Omg4>3.141){
       printf("No, Omega_4 cannot be larger than 3.14. Enter Omega_4 again:\n");
       scanf("%f", &Omg4);
    }
    printf("Enter eta_{-3}: \n");
    scanf("%f", &eta3);
    printf("Enter k_{bow}c_f: \n");
    scanf("%f", &kc);

    printf("\n");
    printf("2. Setting up radio parameters \n");
    printf("Enter L_va/10^29 erg s-1 Hz-1: \n");
    scanf("%f", &Lva);
    if(Lva>1.0e6) {
     printf("L_va should be in units of 10^29 erg s-1 Hz-1. Enter L_va again: \n");
     scanf("%f", &Lva);
    }
    printf("Enter nu_a/GHz : \n");
    scanf("%f", &va);
    if(va>1.0e3) {
     printf("nu_a should be in units of GHz. Enter nu_a again: \n");
     scanf("%f", &va);
    }
    printf("Enter t_2 (t/100days, time since outburst in source rest frame): \n");
    scanf("%f", &t2);
    printf("Enter power-law index p of CRe (positive value. dNe/dgamma~gamma^{-p}): \n");
    scanf("%f", &p);
    if(p<0.0) {
     printf("p should be positive, enter p again: \n");
     scanf("%f", &p);
    }

    printf("\n >>> Check the parameters: \n");
    printf(" Omega_4=%5.2f, eta_{-3}=%5.2f, k_{bow}c_f=%5.2f \n",Omg4,eta3,kc);
    printf(" L_va=%5.2fx10^29, nu_a=%5.2fGHz, t_2=%5.2f, index_p=%5.2f  \n",Lva,va,t2,p);
    vat2=va*t2;

/* Calculating \Gamma listed in Table 1 and presented in Equation 16-19 */
    gm1 =(p+4.0)/(2.0*p+13.0);
    gm2 =-(9.0*p+46.0)/(p+4.0)/(2.0*p+13.0);
    gm3 =-(2.0*p*p+14.0*p+19.0)/(p+4.0)/(2.0*p+13.0);
    gm4 =(p+9.0)/(2.0*p+13.0);
    gm5 =3.0*(p+4.0)/(2.0*p+13.0);
    gm6 =(-11.0*p-34.0)/(p+4.0)/(2.0*p+13.0);
    gm7 =-(2.0*p*p+8.0*p+5.0)/(p+4.0)/(2.0*p+13.0);
    gm8 =(1.0-p)/(2.0*p+13.0);
    gm9 =(6.0-p)/(p+4.0)/(2.0*p+13.0);
    gm10=(3.0*p+7.0)/(p+4.0)/(2.0*p+13.0);

/* Calculating G1,G2,G3 in Equation 19 */
    G0=pow(p-1.5,-2.2)-0.135;
    G1=pow(p/2.5, 1.5)*pow(G0,-gm1);
    G2=pow(p/2.5,-1.5)*pow(G0,-gm5);
    G3=pow(p/2.5,-1.5)*pow(G0,-gm1);

/* Calculating mass outflow rate, kinetic luminosity, and velocity (Equation 16-18)  */
    m  =pow(Lva/1.3,gm1)*vat2/36.0*G1 *pow(kc,-gm1)*pow(eta3,gm2)*pow(Omg4, gm4);
    m  =m  *pow(epsb2,gm3)*pow(epse1,gm2);

    mv2=pow(Lva/1.3,gm5)*36.0/vat2*G2 *pow(kc,-gm5)*pow(eta3,gm6)*pow(Omg4, gm8);
    mv2=mv2*pow(epsb2,gm7)*pow(epse1,gm6);

    v9 =pow(Lva/1.3,gm1)*36.0/vat2*G3 *pow(kc,-gm1)*pow(eta3,gm9)*pow(Omg4,-gm1);
    v9 =v9 *pow(epsb2,gm10)*pow(epse1,gm9);


    Vdn=sqrt(16.0/(1836.0* epse1*0.1)); /* V_DN=sqrt{16m_e/(m_p epbe)} in units of c*/
    printf("\n The critical bulk velocity V_DN=%6.3f c for deep-Newtonian or non-deep-Newtonian \n \n",Vdn);

    printf("======================================================== \n");

/* Now check the result. 1,superluminal or not. 2, V<V_DN or V>V_DN */
    if(v9/30.0 > 1.0) {
      printf("Superluminal! Invalid results. Check the parameters again. \n");
      printf("The current model only applies to non-relativistic outflow.\n");
      return 0;
    }

/* Print the results */
    if(v9/30.0 <= Vdn) {
    printf(">>> Outflow-Cloud Model (deep-Newtonian regime: V < V_DN)  \n");
    printf(">>> The parameters of TDE outflow are: \n \n");
    printf("dotM=%6.3f Msun/yr, Lkin= %6.2fx10^43 erg/s, V=%6.3f c \n \n",m, mv2*3.2,v9/30.0);
    printf("============ Calculated by Equation 16-18. ============= \n");
    }

    if(v9/30.0 > Vdn) {
      printf("\n");
      printf("Caution: V > V_DN. \n");
      printf("Calculation for V>V_DN case is started automatically: \n");
    }

/* Non-deep-Newtonian regime (V>V_DN)  */
/* Calculating \Gm^{'} listed in Table 1 and presented in Equation 20-22 */
    p9=4.0*p+9.0;
    p18=18.0-7.0*p;

    Gms1 =(p+6.0)/p9;
    Gms2 =-(2.0*p+13.0)/p9;
    Gms3 =1.0/p9;
    Gms4 =(-p*p-47.0*p+4.0)/p9/(p+4.0);
    Gms5 =p18/p9;
    Gms6 =2.0*(p+4.0)/(p+6.0)+Gms2*p18/(p+6.0);
    Gms7 =-Gms1*p18/(p+6.0);
    Gms8 =-4.0/(p+6.0)-Gms3*p18/(p+6.0);
    Gms9 =-(p+2.0)/(p+6.0)+Gms3*p18/(p+6.0);
    Gms10=4.0*(1.0-p)/(p+6.0)+Gms4*p18/(p+6.0);
    Gms11=1.0-Gms1*p18/(p+6.0);

/*  G4 and G5 in equation 22 */
    GP4=pow(p/2.5,3.1*Gms1);
    GP5=pow(p/2.5,3.1*Gms5);

/* Calculating velocity and mass outflow rate (Equation 20-21)  */
    beta=pow(Lva/3.4e2,Gms1)*pow(vat2/16.4,Gms2)*GP4*pow(kc,-Gms1)*pow(eta3,-Gms3)*pow(Omg4,-Gms1);
    beta=beta*pow(epsb2,Gms3)*pow(epse1,Gms4);

    dtm =pow(Lva/3.4e2,Gms5)*pow(vat2/16.4,Gms6)*GP5*pow(kc, Gms7)*pow(eta3, Gms8)*pow(Omg4,Gms11);
    dtm =dtm*pow(epsb2,Gms9)*pow(epse1,Gms10);

/* Print the results */
    if(v9/30.0 > Vdn) {
    printf(">>> Outflow-Cloud Model (V > V_DN regime) \n");
    printf(">>> The parameters of TDE outflow are: \n \n");
    printf("dotM=%6.3f Msun/yr, Lkin= %5.2fx10^43 erg/s, Beta=%6.3f \n \n",dtm,2.9e3*dtm*beta*beta,beta);
    printf("============ Calculated by Equation 20,21 ============== \n");
    }

    printf("===========  Written by G.Mou, 2021.Dec.02  ============ \n");
    printf("=============== <<  End of Program  >> ================= \n");

}

