MODULE helper

USE shared_data

CONTAINS


  FUNCTION MomentumFromTemperature(mass,temperature,idum)

    REAL(num), INTENT(IN) :: mass,temperature
    INTEGER, INTENT(INOUT) :: idum
    REAL(num) :: MomentumFromTemperature

    REAL(num) :: mean=0.0_num
    REAL(num) :: stdev
    REAL(num) :: rand1,rand2,w

    !This is a basic polar Box-Muller transform
    !It generates gaussian distributed random numbers
    !The standard deviation (stdev) is related to temperature

    stdev=SQRT(temperature*kb*mass)

    DO
       rand1=Random(idum)
       rand2=Random(idum)

       rand1=2.0_num*rand1 - 1.0_num
       rand2=2.0_num*rand2 - 1.0_num

       w=rand1**2 + rand2**2

       IF (w .LT. 1.0_num) EXIT
    ENDDO

    w = SQRT((-2.0_num * LOG(w) )/w)


    MomentumFromTemperature = rand1 * w * stdev

  END FUNCTION MomentumFromTemperature

  FUNCTION Random(idum)

    INTEGER,INTENT(INOUT) :: idum
    REAL(num) :: Random
    INTEGER,PARAMETER :: IA=16807, IM=2147483647, IQ=127773
    INTEGER,PARAMETER :: IR=2836, MASK=123459876
    REAL(num),PARAMETER :: AM=1.0/REAL(IM,8)

    INTEGER :: k

    idum=XOR(idum,mask)
    k=idum/IQ

    idum=IA*(idum-k*IQ)-IR*k
    IF (idum .LT. 0) idum=idum+IM

    Random=AM*idum
    idum=XOR(idum,mask)

  END FUNCTION Random

END MODULE helper
