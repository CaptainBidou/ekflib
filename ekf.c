/**
 * \file ekf.C
 * \details Code source de ekf.h
 * \author Tomás Salvado Robalo
 * \date 2024
*/

/* ------------------------------------------------------------------------ */
/*                   E N T Ê T E S    S T A N D A R D S                     */
/* ------------------------------------------------------------------------ */
#include "ekf.h"
/* ------------------------------------------------------------------------ */
/*              D É F I N I T I O N S   D E   T Y P E S                     */
/* ------------------------------------------------------------------------ */
/**
 * \struct dg_t
 * \brief Structure représentant un vecteur de dérivées
 * 
 */
typedef struct {
    double mat[3];
}dg_t;


/* ------------------------------------------------------------------------ */
/*            P R O T O T Y P E S    D E    F O N C T I O N S               */
/* ------------------------------------------------------------------------ */
double h(ekf_t * ekf,measure_t * mesure);

double h_derivative(ekf_t * ekf,measure_t * mesure);

double R0(ekf_t * ekf,measure_t * mesure);

double R0_derivative(ekf_t * ekf,measure_t * mesure);

double g(ekf_t * ekf,measure_t * mesure);

dg_t* g_derivative(ekf_t * ekf,measure_t * mesure);

ekf_t * ekfPredict(ekf_t * ekf,measure_t * mesure);

ekf_t * ekfCorrect(ekf_t * ekf,measure_t * mesure);

/* ------------------------------------------------------------------------ */
/*                  C O D E    D E S    F O N C T I O N S                   */
/* ------------------------------------------------------------------------ */

/**
 * \fn void ekf_init(ekf_t * ekf)
 * \brief Initialisation de l'ekf
 * \param ekf_t * ekf pointeur sur l'ekf
 */
ekf_t * ekf_init(void){
    //allocation de la mémoire pour l'ekf
    ekf_t * ekf = (ekf_t *)malloc(sizeof(ekf_t));

    //initialisation de la matrice d'état
    ekf->x[0][0] = 0;
    ekf->x[1][0] = 0;
    ekf->x[2][0] = 0.5;

    //initialisation de la matrice de covariance de l'état
    ekf->P[0][0] = 0.001;
    ekf->P[0][1] = 0;
    ekf->P[0][2] = 0;
    ekf->P[1][0] = 0;
    ekf->P[1][1] = 0.001;
    ekf->P[1][2] = 0;
    ekf->P[2][0] = 0;
    ekf->P[2][1] = 0;
    ekf->P[2][2] = 0.01;

    //initialisation de la matrice de covariance du bruit de processus
    ekf->Q[0][0] = 0.000001;
    ekf->Q[0][1] = 0;
    ekf->Q[0][2] = 0;
    ekf->Q[1][0] = 0;
    ekf->Q[1][1] = 0.000001;
    ekf->Q[1][2] = 0;
    ekf->Q[2][0] = 0;
    ekf->Q[2][1] = 0;
    ekf->Q[2][2] = 0.0001;

    //initialisation de la matrice de covariance du bruit de mesure
    ekf->R = 0.2;

    //initialisation de la matrice de transition de l'état
    ekf->A[0][0] = 1-(SAMPLE_RATE/(R1*C1));
    ekf->A[0][1] = 0;
    ekf->A[0][2] = 0;
    ekf->A[1][0] = 0;
    ekf->A[1][1] = 1-(SAMPLE_RATE/(R2*C2));
    ekf->A[1][2] = 0;
    ekf->A[2][0] = 0;
    ekf->A[2][1] = 0;
    ekf->A[2][2] = 1;

    //initialisation de la matrice de contrôle de l'état
    ekf->B[0][0] = SAMPLE_RATE/C1;
    ekf->B[1][0] = SAMPLE_RATE/C2;
    ekf->B[2][0] = -SAMPLE_RATE/(3600*QN);

    ekf->K[0][0] = 0;
    ekf->K[1][0] = 0;
    ekf->K[2][0] = 0;

    return ekf;
}

/**
 * \fn void ekf_run(ekf_t * ekf)
 * \brief Exécution de l'ekf
 * \param ekf_t * ekf pointeur sur l'ekf
 * \param measure_t * mesure pointeur sur la mesure
 */
estimation_t* ekf_run(ekf_t * ekf,measure_t * mesure){
    //prédiction
    ekf = ekfPredict(ekf,mesure);
    //correction
    ekf = ekfCorrect(ekf,mesure);

    //allocation de la mémoire pour l'estimation
    estimation_t * estimation = (estimation_t *)malloc(sizeof(estimation_t));

    //estimation de l'état de charge
    estimation->soc = ekf->x[2][0];
    estimation->voltage = g(ekf,mesure);

    return estimation;

}

/**
 * \fn double h(ekf_t * ekf,measure_t * mesure)
 * \brief fonction de mesure
 * \param ekf_t * ekf pointeur sur l'ekf
 * \param measure_t * mesure pointeur sur la mesure
 */
double h(ekf_t * ekf,measure_t * mesure){
    double x = ekf->x[2][0];
    return pow(10,4)*(pow(x,9)*0.140437566214432+\
    pow(x,8)*-0.693293580320924+pow(x,7)*1.448317451181394 + \
    pow(x,6) *-1.665299094951629 + pow(x,5)*1.148704101226141 + \
    pow(x,4)*-0.486836353839831 + pow(x,3)*0.125420712206318 + \
    pow(x,2)*-0.018961803736654 + x*0.001657801378501 + 0.000269333059573);
}

/**
 * \fn double h_derivative(ekf_t * ekf,measure_t * mesure)
 * \brief dérivée de la fonction de mesure
 * \param ekf_t * ekf pointeur sur l'ekf
 * \param measure_t * mesure pointeur sur la mesure
 */
double h_derivative(ekf_t * ekf,measure_t * mesure){
    double x = ekf->x[2][0];
    return pow(10,4)*(9*pow(x,8)*0.140437566214432+\
    8*pow(x,7)*-0.693293580320924+7*pow(x,6)*1.448317451181394 + \
    6*pow(x,5)*-1.665299094951629 + 5*pow(x,4)*1.148704101226141 + \
    4*pow(x,3)*-0.486836353839831 + 3*pow(x,2)*0.125420712206318 + \
    2*pow(x,1)*-0.018961803736654 + 1*0.001657801378501);
}


/**
 * \fn double R0(ekf_t * ekf,measure_t * mesure)
 * \brief fonction de mesure
 * \param ekf_t * ekf pointeur sur l'ekf
 * \param measure_t * mesure pointeur sur la mesure
 */
double R0(ekf_t * ekf,measure_t * mesure){
    double x = ekf->x[2][0];
    return pow(10,3)*(pow(x,9)*0.440568380336331+\
    pow(x,8)*-2.188575118770938+pow(x,7)*4.662025929324535+\
    pow(x,6)*-5.561277160719505+pow(x,5)*4.069003040512039+\
    pow(x,4)*-1.878727644202677+pow(x,3)*0.541295950462937+\
    pow(x,2)*-0.092097275963785+x*0.008056926448651-0.000160671690337);
}

/**
 * \fn double R0_derivative(ekf_t * ekf,measure_t * mesure)
 * \brief dérivée de la fonction de mesure
 * \param ekf_t * ekf pointeur sur l'ekf
 * \param measure_t * mesure pointeur sur la mesure
 */
double R0_derivative(ekf_t * ekf,measure_t * mesure){
    double x = ekf->x[2][0];
    return pow(10,3)*(9*pow(x,8)*0.440568380336331+\
    8*pow(x,7)*-2.188575118770938+7*pow(x,6)*4.662025929324535+\
    6*pow(x,5)*-5.561277160719505+5*pow(x,4)*4.069003040512039+\
    4*pow(x,3)*-1.878727644202677+3*pow(x,2)*0.541295950462937+\
    2*pow(x,1)*-0.092097275963785+1*0.008056926448651);
}

/**
 * \fn double g(ekf_t * ekf,measure_t * mesure)
 * \brief fonction de mesure
 * \param ekf_t * ekf pointeur sur l'ekf
 * \param measure_t * mesure pointeur sur la mesure
 */
double g(ekf_t * ekf,measure_t * mesure){
    return ekf->x[0][0] - ekf->x[1][0] - R0(ekf,mesure)*mesure->current + h(ekf,mesure);
}

/**
 * \fn double g_derivative(ekf_t * ekf,measure_t * mesure)
 * \brief dérivée de la fonction de mesure
 * \param ekf_t * ekf pointeur sur l'ekf
 * \param measure_t * mesure pointeur sur la mesure
 */
dg_t* g_derivative(ekf_t * ekf,measure_t * mesure){
    dg_t * dg = (dg_t *)malloc(sizeof(dg_t));
    dg->mat[0] = -1;
    dg->mat[1] = -1;
    dg->mat[2] = h_derivative(ekf,mesure) - R0_derivative(ekf,mesure)*mesure->current;
    return dg;

}

/**
 * \fn ekf_t * ekfPredict(ekf_t * ekf,measure_t * mesure)
 * \brief prédiction de l'ekf
 * \param ekf_t * ekf pointeur sur l'ekf
 * \param measure_t * mesure pointeur sur la mesure
 */
ekf_t * ekfPredict(ekf_t * ekf,measure_t * mesure){
    //prédiction de l'état
    ekf->x[0][0] = ekf->x[0][0]*(ekf->A[0][0]+ekf->A[0][1]+ekf->A[0][2]) + ekf->B[0][0]*mesure->current;
    ekf->x[1][0] = ekf->x[1][0]*(ekf->A[1][0]+ekf->A[1][1]+ekf->A[1][2]) + ekf->B[1][0]*mesure->current;
    ekf->x[2][0] = ekf->x[2][0]*(ekf->A[2][0]+ekf->A[2][1]+ekf->A[2][2]) + ekf->B[2][0]*mesure->current;
    //prédiction de la covariance de l'état 
    double P[3][3];
    P[0][0] = ekf->A[0][0]*ekf->P[0][0]+ekf->A[0][1]*ekf->P[1][0]+ekf->A[0][2]*ekf->P[2][0];
    P[0][1] = ekf->A[0][0]*ekf->P[0][1]+ekf->A[0][1]*ekf->P[1][1]+ekf->A[0][2]*ekf->P[2][1];
    P[0][2] = ekf->A[0][0]*ekf->P[0][2]+ekf->A[0][1]*ekf->P[1][2]+ekf->A[0][2]*ekf->P[2][2];
    P[1][0] = ekf->A[1][0]*ekf->P[0][0]+ekf->A[1][1]*ekf->P[1][0]+ekf->A[1][2]*ekf->P[2][0];
    P[1][1] = ekf->A[1][0]*ekf->P[0][1]+ekf->A[1][1]*ekf->P[1][1]+ekf->A[1][2]*ekf->P[2][1];
    P[1][2] = ekf->A[1][0]*ekf->P[0][2]+ekf->A[1][1]*ekf->P[1][2]+ekf->A[1][2]*ekf->P[2][2];
    P[2][0] = ekf->A[2][0]*ekf->P[0][0]+ekf->A[2][1]*ekf->P[1][0]+ekf->A[2][2]*ekf->P[2][0];
    P[2][1] = ekf->A[2][0]*ekf->P[0][1]+ekf->A[2][1]*ekf->P[1][1]+ekf->A[2][2]*ekf->P[2][1];
    P[2][2] = ekf->A[2][0]*ekf->P[0][2]+ekf->A[2][1]*ekf->P[1][2]+ekf->A[2][2]*ekf->P[2][2];

    ekf->P[0][0] = P[0][0]*ekf->A[0][0]+P[0][1]*ekf->A[0][1]+P[0][2]*ekf->A[0][2]+ekf->Q[0][0];
    ekf->P[0][1] = P[0][0]*ekf->A[0][0]+P[0][1]*ekf->A[0][1]+P[0][2]*ekf->A[0][2]+ekf->Q[0][1];
    ekf->P[0][2] = P[0][0]*ekf->A[0][0]+P[0][1]*ekf->A[0][1]+P[0][2]*ekf->A[0][2]+ekf->Q[0][2];
    ekf->P[1][0] = P[1][0]*ekf->A[1][0]+P[1][1]*ekf->A[1][1]+P[1][2]*ekf->A[1][2]+ekf->Q[1][0];
    ekf->P[1][1] = P[1][0]*ekf->A[1][0]+P[1][1]*ekf->A[1][1]+P[1][2]*ekf->A[1][2]+ekf->Q[1][1];
    ekf->P[1][2] = P[1][0]*ekf->A[1][0]+P[1][1]*ekf->A[1][1]+P[1][2]*ekf->A[1][2]+ekf->Q[1][2];
    ekf->P[2][0] = P[2][0]*ekf->A[2][0]+P[2][1]*ekf->A[2][1]+P[2][2]*ekf->A[2][2]+ekf->Q[2][0];
    ekf->P[2][1] = P[2][0]*ekf->A[2][0]+P[2][1]*ekf->A[2][1]+P[2][2]*ekf->A[2][2]+ekf->Q[2][1];
    ekf->P[2][2] = P[2][0]*ekf->A[2][0]+P[2][1]*ekf->A[2][1]+P[2][2]*ekf->A[2][2]+ekf->Q[2][2];
    
    return ekf;
}

/**
 * \fn ekf_t * ekfCorrect(ekf_t * ekf,measure_t * mesure)
 * \brief correction de l'ekf
 * \param ekf_t * ekf pointeur sur l'ekf
 * \param measure_t * mesure pointeur sur la mesure
 */
ekf_t * ekfCorrect(ekf_t * ekf,measure_t * mesure){
    
    dg_t* dg = g_derivative(ekf,mesure);
    double index0 = dg->mat[0]*ekf->P[0][0]+dg->mat[1]*ekf->P[1][0]+dg->mat[2]*ekf->P[2][0];
    double index1 = dg->mat[0]*ekf->P[0][1]+dg->mat[1]*ekf->P[1][1]+dg->mat[2]*ekf->P[2][1];
    double index2 = dg->mat[0]*ekf->P[0][2]+dg->mat[1]*ekf->P[1][2]+dg->mat[2]*ekf->P[2][2];

    double gain = index0*dg->mat[0]+index1*dg->mat[1]+index2*dg->mat[2]+ekf->R;
    

    if(gain > DBL_MAX/1000){
        printf("gain = 0\n");
        gain = 0;
    }else{
        gain = 1/gain;
        printf("gain = %f\n",gain);
    }

    //calcul du gain de Kalman
    ekf->K[0][0] = (ekf->P[0][0]*dg->mat[0]+ekf->P[0][1]*dg->mat[1]+ekf->P[0][2]*dg->mat[2])*gain;
    ekf->K[1][0] = (ekf->P[1][0]*dg->mat[0]+ekf->P[1][1]*dg->mat[1]+ekf->P[1][2]*dg->mat[2])*gain;
    ekf->K[2][0] = (ekf->P[2][0]*dg->mat[0]+ekf->P[2][1]*dg->mat[1]+ekf->P[2][2]*dg->mat[2])*gain;

    double gResult = g(ekf,mesure);

    //correction de l'état
    ekf->x[0][0] = ekf->x[0][0] + ekf->K[0][0]*(mesure->voltage - gResult);
    ekf->x[1][0] = ekf->x[1][0] + ekf->K[1][0]*(mesure->voltage - gResult);
    ekf->x[2][0] = ekf->x[2][0] + ekf->K[2][0]*(mesure->voltage - gResult);

    //correction de la covariance de l'état
    dg = g_derivative(ekf,mesure);
    double ikdg[3][3] = {1-ekf->K[0][0]*dg->mat[0],-ekf->K[0][0]*dg->mat[1],-ekf->K[0][0]*dg->mat[2],
                         -ekf->K[1][0]*dg->mat[0],1-ekf->K[1][0]*dg->mat[1],-ekf->K[1][0]*dg->mat[2],
                         -ekf->K[2][0]*dg->mat[0],-ekf->K[2][0]*dg->mat[1],1-ekf->K[2][0]*dg->mat[2]};

    double P[3][3];

    P[0][0] = ikdg[0][0]*ekf->P[0][0]+ikdg[0][1]*ekf->P[1][0]+ikdg[0][2]*ekf->P[2][0];
    P[0][1] = ikdg[0][0]*ekf->P[0][1]+ikdg[0][1]*ekf->P[1][1]+ikdg[0][2]*ekf->P[2][1];
    P[0][2] = ikdg[0][0]*ekf->P[0][2]+ikdg[0][1]*ekf->P[1][2]+ikdg[0][2]*ekf->P[2][2];
    P[1][0] = ikdg[1][0]*ekf->P[0][0]+ikdg[1][1]*ekf->P[1][0]+ikdg[1][2]*ekf->P[2][0];
    P[1][1] = ikdg[1][0]*ekf->P[0][1]+ikdg[1][1]*ekf->P[1][1]+ikdg[1][2]*ekf->P[2][1];
    P[1][2] = ikdg[1][0]*ekf->P[0][2]+ikdg[1][1]*ekf->P[1][2]+ikdg[1][2]*ekf->P[2][2];
    P[2][0] = ikdg[2][0]*ekf->P[0][0]+ikdg[2][1]*ekf->P[1][0]+ikdg[2][2]*ekf->P[2][0];
    P[2][1] = ikdg[2][0]*ekf->P[0][1]+ikdg[2][1]*ekf->P[1][1]+ikdg[2][2]*ekf->P[2][1];
    P[2][2] = ikdg[2][0]*ekf->P[0][2]+ikdg[2][1]*ekf->P[1][2]+ikdg[2][2]*ekf->P[2][2];

    ekf->P[0][0] = P[0][0];
    ekf->P[0][1] = P[0][1];
    ekf->P[0][2] = P[0][2];
    ekf->P[1][0] = P[1][0];
    ekf->P[1][1] = P[1][1];
    ekf->P[1][2] = P[1][2];
    ekf->P[2][0] = P[2][0];
    ekf->P[2][1] = P[2][1];
    ekf->P[2][2] = P[2][2];

    P[0][0] = ekf->P[0][0]*ikdg[0][0]+ekf->P[0][1]*ikdg[0][1]+ekf->P[0][2]*ikdg[0][2] + ekf->K[0][0]*ekf->K[0][0]*ekf->R;
    P[0][1] = ekf->P[0][0]*ikdg[0][0]+ekf->P[0][1]*ikdg[0][1]+ekf->P[0][2]*ikdg[0][2] + ekf->K[0][0]*ekf->K[1][0]*ekf->R;
    P[0][2] = ekf->P[0][0]*ikdg[0][0]+ekf->P[0][1]*ikdg[0][1]+ekf->P[0][2]*ikdg[0][2] + ekf->K[0][0]*ekf->K[2][0]*ekf->R;
    P[1][0] = ekf->P[1][0]*ikdg[1][0]+ekf->P[1][1]*ikdg[1][1]+ekf->P[1][2]*ikdg[1][2] + ekf->K[1][0]*ekf->K[0][0]*ekf->R;
    P[1][1] = ekf->P[1][0]*ikdg[1][0]+ekf->P[1][1]*ikdg[1][1]+ekf->P[1][2]*ikdg[1][2] + ekf->K[1][0]*ekf->K[1][0]*ekf->R;
    P[1][2] = ekf->P[1][0]*ikdg[1][0]+ekf->P[1][1]*ikdg[1][1]+ekf->P[1][2]*ikdg[1][2] + ekf->K[1][0]*ekf->K[2][0]*ekf->R;
    P[2][0] = ekf->P[2][0]*ikdg[2][0]+ekf->P[2][1]*ikdg[2][1]+ekf->P[2][2]*ikdg[2][2] + ekf->K[2][0]*ekf->K[0][0]*ekf->R;
    P[2][1] = ekf->P[2][0]*ikdg[2][0]+ekf->P[2][1]*ikdg[2][1]+ekf->P[2][2]*ikdg[2][2] + ekf->K[2][0]*ekf->K[1][0]*ekf->R;
    P[2][2] = ekf->P[2][0]*ikdg[2][0]+ekf->P[2][1]*ikdg[2][1]+ekf->P[2][2]*ikdg[2][2] + ekf->K[2][0]*ekf->K[2][0]*ekf->R;

    ekf->P[0][0] = P[0][0];
    ekf->P[0][1] = P[0][1];
    ekf->P[0][2] = P[0][2];
    ekf->P[1][0] = P[1][0];
    ekf->P[1][1] = P[1][1];
    ekf->P[1][2] = P[1][2];
    ekf->P[2][0] = P[2][0];
    ekf->P[2][1] = P[2][1];
    ekf->P[2][2] = P[2][2];




    return ekf;

}