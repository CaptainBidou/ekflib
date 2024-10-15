/**
 * \file ekf.h
 * \details déclaration des fonctions et types de ekf.c
 * \author Tomás Salvado Robalo
 * \date 2024
*/
#ifndef EKF_H
#define EKF_H
#endif

/* ------------------------------------------------------------------------ */
/*                   E N T Ê T E S    S T A N D A R D S                     */
/* ------------------------------------------------------------------------ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <limits.h>



/* ------------------------------------------------------------------------ */
/*              C O N S T A N T E S     S Y M B O L I Q U E S               */
/* ------------------------------------------------------------------------ */
#define QN 3.08 //Valeur de la charge de la cellule
#define SAMPLE_RATE 1 //Fréquence d'échantillonage
#define R1 0.009818410271493//Valeur de la résistance de la cellule
#define C1 15639.54740330107
#define R2 0.031463057438081
#define C2 3933.292323912280

/* ------------------------------------------------------------------------ */
/*                    M A C R O    F O N C T I O N S                        */
/* ------------------------------------------------------------------------ */



/* ------------------------------------------------------------------------ */
/*              D É F I N I T I O N S   D E   T Y P E S                     */
/* ------------------------------------------------------------------------ */

/**
 * \struct ekf_t
 * \brief Structure représentant un filtre de Kalman étendu
 * 
 * La structure contient les matrices et vecteurs nécessaires pour le filtre de Kalman étendu.
*/
typedef struct {
    double x[3][1]; /*!< Vecteur d'état */
    double P[3][3]; /*!< Matrice de covariance de l'état */
    double Q[3][3]; /*!< Matrice de covariance du bruit de processus */
    double R; /*!< Matrice de covariance du bruit de mesure */
    double A[3][3]; /*!< Matrice de transition de l'état */
    double B[3][1]; /*!< Matrice de contrôle de l'état */
    double K[3][1]; /*!< Matrice de gain de Kalman */
} ekf_t;

/**
 * \struct estimation_t
 * \brief Structure représentant une estimation
 * 
 * La structure contient les valeurs de l'état de charge et de la tension.
*/
typedef struct {
    double soc; /*!< Valeur de l'état de charge */
    double voltage; /*!< valeur de la tension */
}estimation_t;

/**
 * \struct mesure_t
 * \brief Structure représentant une mesure
 * 
 * La structure contient les valeurs de la tension et du courant.
 * Lorsque le courant est négatif cela signifie que la batterie se charge 
 * et lorsqu'il est positif cela signifie que la batterie se décharge.
*/
typedef struct{
    double voltage; /*!< Valeur de la tension */
    double current; /*!< Valeur du courant */
}mesure_t;

/* ------------------------------------------------------------------------ */
/*            P R O T O T Y P E S    D E    F O N C T I O N S               */
/* ------------------------------------------------------------------------ */

/**
 * \fn ekf_t * ekf_init(void)
 * \brief Initialisation de l'ekf
 * \return ekf_t * pointeur sur l'ekf
 */
ekf_t * ekf_init(void);

/**
 * \fn estimation_t* ekf_run(ekf_t * ekf,mesure_t * mesure)
 * \brief Exécution de l'ekf
 * \param ekf_t * ekf pointeur sur l'ekf
 * \param mesure_t * mesure pointeur sur la mesure
 * \return estimation_t* pointeur sur l'estimation
 */
estimation_t* ekf_run(ekf_t * ekf,mesure_t * mesure);