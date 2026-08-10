#imports: Numpy for mathematical operations.
import numpy as np

#Calculates the time derivatives of spin amplitudes using Schrödingers equation.
#Inputs: Alpha and beta are spin amplitudes. Bx, By, Bz are magnetic field components in each direction.
#Outputs: Dalpha and dbeta are the time derivatives.
def aBDerivative(alpha, beta, Bx, By, Bz):
    dalpha = (((-1.0j)*Bz)*alpha) + ((((-1.0j)*Bx)-(By))*beta)
    dbeta = ((((-1.0j)*Bx)+(By))*alpha) + (((1.0j)*Bz)*beta)
    return dalpha, dbeta

#Calculates the spin state for the next time step using Euler's method.
#Inputs: Alpha and beta are spin amplitudes. Bx, By, Bz are magnetic field components in each direction. Deltat is the timestep.
#Outputs: Outputs the updated spin amplitudes at time t + deltat.
def euler(alpha, beta, Bx, By, Bz, deltat):
    dalpha, dbeta = aBDerivative(alpha, beta, Bx, By, Bz)
    return ((dalpha*deltat)+alpha), ((dbeta*deltat)+beta)

#Calculates the spin state for the next time step using Improved Euler's method or Heun method.
#Inputs: Alpha and beta are spin amplitudes. Bx, By, Bz are magnetic field components in each direction. Deltat is the timestep.
#Outputs: Outputs the updated spin amplitudes at time t + deltat.
def improvedEuler(alpha, beta, Bx, By, Bz, deltat):
    alphak1, betak1 = aBDerivative(alpha, beta, Bx, By, Bz)
    alphag = alpha + (alphak1*deltat)
    betag = beta + (betak1*deltat)
    alphak2, betak2 = aBDerivative(alphag, betag, Bx, By, Bz)
    return (alpha + (0.5*deltat)*(alphak1 + alphak2)), (beta + (0.5*deltat)*(betak1 + betak2))

#Calculates the spin state for the next time step using the Fourth Order Runge-Kutta method.
#Inputs: Alpha and beta are spin amplitudes. Bx and By are magnetic field components, BzList contains the time-dependent Bz values, deltat is the timestep, and i is the current time index.
#Outputs: Outputs the updated spin amplitudes at time t + deltat.
def rk4(alpha, beta, Bx, By, BzList, deltat, i):
    alphak1, betak1 = aBDerivative(alpha, beta, Bx, By, BzList[i])
    midalpha12 = alpha+((deltat/2)*(alphak1))
    midbeta12 = beta+((deltat/2)*(betak1))
    alphak2, betak2 = aBDerivative(midalpha12, midbeta12, Bx, By, (((BzList[i+1]+BzList[i]))*0.5))
    midalpha23 = alpha+((deltat/2)*(alphak2))
    midbeta23 = beta+((deltat/2)*(betak2))
    alphak3, betak3 = aBDerivative(midalpha23, midbeta23, Bx, By, (((BzList[i+1]+BzList[i]))*0.5))
    midalpha34 = alpha+((deltat)*(alphak3))
    midbeta34 = beta+((deltat)*(betak3))
    alphak4, betak4 = aBDerivative(midalpha34, midbeta34, Bx, By, BzList[i+1])
    return (alpha+((deltat/6)*(alphak1+(2*alphak2)+(2*alphak3)+alphak4))), (beta+((deltat/6)*(betak1+(2*betak2)+(2*betak3)+betak4)))

#Calculates the exact spin state after evolving under a constant magnetic field.
#Inputs: Alpha and beta are spin amplitudes. Bx, By, Bz are constant magnetic field components. t is the evolution time.
#Outputs: Outputs the spin amplitudes after evolving for time t.
def constantBSolution(alpha, beta, Bx, By, Bz, t):
    B=np.sqrt((Bx*Bx)+(By*By)+(Bz*Bz))
    if(B == 0):
        return alpha, beta
    cosPart = np.array([[np.cos(B*t),0],[0,np.cos(B*t)]])
    sinPart = np.array([[Bz, (Bx-((1.0j)*(By)))], [Bx+((1.0j)*By), -Bz]])
    sinPart = sinPart * ((-1.0j)*(1/B)*(np.sin(B*t)))
    result = sinPart + cosPart
    initialArray = np.array([alpha, beta])
    finalArray = result @ initialArray
    return (finalArray[0]), (finalArray[1])

#Evolves a singular spin trajectory over the entire simulation based on which numerical method was selected.
#Inputs: Alpha and beta are the initial amplitudes, Bx and By are the components of the magnetic field, BzList is the time-dependent magnetic field, deltat is the time step, tmax is the full time of the evolution, and method is the chosen numerical method. 
#Outputs: List of every time the system is evolved over, all alpha values, all beta values, all expectation values of the x-component of the spin, and the normalization error.
def solution(alpha, beta, Bx, By, BzList, deltat, tmax, method):
    timeList = []
    alphaList = []
    betaList = []
    xList = []
    normList = []
    for i in range((int(tmax/deltat))+1):
        timeList.append(deltat*i)
        alphaList.append(alpha)
        betaList.append(beta)
        xList.append(2*(((alpha.conjugate())*(beta)).real))
        normList.append((((abs(alpha))*(abs(alpha))+(abs(beta))*(abs(beta)))-1))
        if(i < (int(tmax/deltat))):
            if(method == 1):
                alpha, beta = euler(alpha, beta, Bx, By, BzList[i], deltat)
            elif(method == 2):
                alpha, beta = improvedEuler(alpha, beta, Bx, By, BzList[i], deltat)
            elif(method == 3):
                alpha, beta = rk4(alpha, beta, Bx, By, BzList, deltat, i)
            elif(method == 4):
                alpha, beta = constantBSolution(alpha, beta, Bx, By, BzList[i], deltat)
    return timeList, alphaList, betaList, xList, normList

#Simulates lots of Bz trajectories and averages the spin evolution.
#Inputs: Alpha and beta are the initial amplitudes, Bx and By are magnetic field components, BzData contains multiple Bz trajectories, deltat is the time step, tmax is the full evolution time, and method is the chosen numerical method.
#Outputs: List of Bz trajectories, time values, and the averaged alpha, beta, X expectation values, and normalization error.
def Bzaverages(alpha, beta, Bx, By, BzData, deltat, tmax, method):
    BzTrajectoryList = []
    newTimeList = None
    newAlphaList = None
    newBetaList = None
    newXList = None
    newNormList = None
    for i in range((len(BzData))):
        BzList = BzData[i]
        timeList, alphaList, betaList, xList, normList = solution(alpha, beta, Bx, By, BzList, deltat, tmax, method)
        BzTrajectoryList.append(BzData[i])
        if(i == 0):
            newTimeList = timeList
            newAlphaList = np.array(alphaList)
            newBetaList = np.array(betaList)
            newXList = np.array(xList)
            newNormList = np.array(normList)
        else:
            newAlphaList = newAlphaList+np.array(alphaList)
            newBetaList = newBetaList+np.array(betaList)
            newXList = newXList+np.array(xList)
            newNormList = newNormList+np.array(normList)
    newAlphaList /= (len(BzData))
    newBetaList /= (len(BzData))
    newXList /= (len(BzData))
    newNormList /= (len(BzData))
    return BzTrajectoryList, newTimeList, newAlphaList, newBetaList, newXList, newNormList