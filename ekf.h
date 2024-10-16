/**
 * \file ekf.h
 * \details Declaration of functions and types from ekf.c
 * \author Tomás Salvado Robalo
 * \date 2024
*/
#ifndef EKF_H
#define EKF_H


/* ------------------------------------------------------------------------ */
/*                      S T A N D A R D    H E A D E R S                    */
/* ------------------------------------------------------------------------ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <float.h>
#include <limits.h>



/* ------------------------------------------------------------------------ */
/*                    S Y M B O L I C    C O N S T A N T S                  */
/* ------------------------------------------------------------------------ */
#define QN 3.08 //Cell charge value
#define SAMPLE_RATE 1 //Sampling frequency
#define R1 0.009818410271493 //Cell resistance value
#define C1 15639.54740330107
#define R2 0.031463057438081
#define C2 3933.292323912280

/* ------------------------------------------------------------------------ */
/*                        M A C R O    F U N C T I O N S                    */
/* ------------------------------------------------------------------------ */



/* ------------------------------------------------------------------------ */
/*                      T Y P E   D E F I N I T I O N S                     */
/* ------------------------------------------------------------------------ */

/**
 * \struct ekf_t
 * \brief Structure representing an extended Kalman filter
 * 
 * The structure contains the matrices and vectors needed for the extended Kalman filter.
*/
typedef struct {
    double x[3][1]; /*!< State vector */
    double P[3][3]; /*!< State covariance matrix */
    double Q[3][3]; /*!< Process noise covariance matrix */
    double R; /*!< Measurement noise covariance matrix */
    double A[3][3]; /*!< State transition matrix */
    double B[3][1]; /*!< State control matrix */
    double K[3][1]; /*!< Kalman gain matrix */
} ekf_t;

/**
 * \struct estimation_t
 * \brief Structure representing an estimation
 * 
 * The structure contains the state of charge and voltage values.
*/
typedef struct {
    double soc; /*!< State of charge value */
    double voltage; /*!< Voltage value */
} estimation_t;

/**
 * \struct mesure_t
 * \brief Structure representing a measurement
 * 
 * The structure contains the voltage and current values.
 * When the current is negative, it indicates the battery is charging,
 * and when it is positive, it indicates the battery is discharging.
*/
typedef struct {
    double voltage; /*!< Voltage value */
    double current; /*!< Current value */
} measure_t;

/* ------------------------------------------------------------------------ */
/*                        F U N C T I O N   P R O T O T Y P E S             */
/* ------------------------------------------------------------------------ */

/**
 * \fn ekf_t * ekf_init(void)
 * \brief Initialization of the ekf
 * \return ekf_t * pointer to the ekf
 */
ekf_t * ekf_init(void);

/**
 * \fn estimation_t* ekf_run(ekf_t * ekf, mesure_t * mesure)
 * \brief Execution of the ekf
 * \param ekf_t * ekf pointer to the ekf
 * \param measure_t * mesure pointer to the measurement
 * \return estimation_t* pointer to the estimation
 */
estimation_t* ekf_run(ekf_t * ekf, measure_t * measure);

#endif
