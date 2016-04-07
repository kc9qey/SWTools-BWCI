!-----------------------------------------------------------------------
!    Module:       numerical_module
!    Author:       S. W. Bruenn
!    Date:         8/16/02
!
!    Common numerical constants
!-----------------------------------------------------------------------

module numerical_module

USE kind_module
SAVE


REAL(KIND=double), PARAMETER  :: zero      = 0.0d0
REAL(KIND=double), PARAMETER  :: half      = 0.5d0
REAL(KIND=double), PARAMETER  :: one       = 1.0d0
REAL(KIND=double), PARAMETER  :: epsilon   = 1.0d-100

END module numerical_module
