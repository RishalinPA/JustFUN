program wh3d
    !use io
    implicit none

	real(8), dimension(:,:), allocatable :: Ne, nu, omega_p2, sigma_0
	real(8) ::  omega_ce
	real(8) :: dBphi_dz, dBr_dz,dBz_dr, drBphi_dr, dEphi_dr
	real(8) :: dEphi_dz, dEr_dz,dEz_dr, drEphi_dr, dBphi_dr
	real(8), dimension(:,:), allocatable :: Er, Ephi, Ez, Br, Bphi, Bz, jext
    real(8), dimension(:,:), allocatable :: jr, jphi, jz
	real(8), dimension(:), allocatable :: View
	real(8), parameter :: pi = 2d0*acos(0d0), f = 158.0d6, c = 3.0d10, Ra = 7.0d0, B0z = 90.0d0
	real(8), parameter :: dr = 0.1d0, dz = 0.2d0, e=4.8d-10, me=.91d-27
	real(8) :: dt, omega, r_i, tau
	real(8) :: startTime, stopTime
	integer(4) :: i, ir, iz, it, r_source, z_source
	integer(4), parameter :: nr = 400, nz = 800, nt = 1000000 !3000000 !расчет на 40 часов
	real(8), parameter :: j = (0.0d0, 1.0d0)
	integer(4) :: int1, int2
	
	!call cpu_time(startTime)
	call system_clock(int1)
	r_source = 20!nint(real(nr-1)/2)
	z_source = 250!nint(real(nz-1)/2)
	print *, 'source:', r_source*dr, z_source*dz!, 'and', r_source*dr+Ra, z_source*dz
	
	!Из Куранта. В цилиндрических без зависимости от фи
	
	dt = 0.05d0 / c / (1.0d0 / dr**2 + 1.0d0 / dz**2)**0.5
	print *, 't_max = ', nt*dt, 'dt=', dt
	!Сделать внешний ввод
	allocate(Er(0:nr-1, 0:nz-1), Ephi(0:nr-1, 0:nz-1), Ez(0:nr-1, 0:nz-1), jext(0:nr-1, 0:nz-1))
	allocate(Br(0:nr-1, 0:nz-1), Bphi(0:nr-1, 0:nz-1), Bz(0:nr-1, 0:nz-1))
	allocate(Ne(0:nr-1, 0:nz-1), nu(0:nr-1, 0:nz-1))
	allocate(jr(0:nr-1, 0:nz-1), jphi(0:nr-1, 0:nz-1), jz(0:nr-1, 0:nz-1), sigma_0(0:nr-1, 0:nz-1))
	allocate(omega_p2(0:nr-1, 0:nz-1), View(nt))
	
	Ne = 3.0d11!reshape([real(8):: (0d0, i = 1, nr*nz)], [nr, nz]) !3d11
	call concentration(Ne, Ne)
	nu = 1.0d6!1.0d6
	
	jz = 0.0d0
	jr = 0.0d0
	jphi = 0.0d0
	
	Er = 0.0d0
	Ephi = 0.0d0
	Ez = 0.0d0
	jext = 0.0d0

	Br = 0.0d0
	Bphi = 0.0d0
	Bz = 0.0d0
	
	dBphi_dz = 0.0d0
	dBr_dz= 0.0d0
	dBz_dr= 0.0d0
	drBphi_dr = 0.0d0
	dBphi_dr = 0.0d0
	dEphi_dz= 0.0d0
	dEr_dz= 0.0d0
	dEphi_dr = 0.0d0
	dEz_dr= 0.0d0
	drEphi_dr = 0.0d0
	sigma_0 = 0.0d0

	View = 0.0d0

	call absorb_layer(f, sigma_0) !поглощающий слой проработать
	call omega_psq(Ne, omega_p2)
	call omega_c(B0z, omega_ce)

	open(23, file = 'sigma.bin', action = 'write', form = 'unformatted')
		write(23) sigma_0
	close(23)

	open(25, file = 'Ne.bin', action = 'write', form = 'unformatted')
		write(25) Ne
	close(25)

	!call exit(0)
    !Пока просто формулы, нет набега фазы
	do it = 0, nt
		!Bz(0, :) = Bz(0, :) - 4* dt * c * Ephi(1, :)/dr !регуляризация в нуле (см Chen1996)

	call source(f, dt, dr, dz, Ra, it, r_source, z_source, jext)

	!$omp parallel workshare
	Ephi = Ephi + jext
	!$omp end parallel workshare

	!$omp parallel do private(dBphi_dr, drBphi_dr, dEphi_dz, dEphi_dr, drEphi_dr, iz)
        do iz=1, nz-1
			do ir = 0, nr-2
				dBphi_dr = (Bphi(ir+1, iz) - Bphi(ir, iz) ) / dr
				drBphi_dr = dBphi_dr + (Bphi(ir+1, iz)  + Bphi(ir, iz))/(2.0d0*(real(ir, 8)+0.5d0) *dr)

				dEphi_dz = (Ephi(ir, iz+1) - Ephi(ir, iz)) / dz
				dEphi_dr = (Ephi(ir+1, iz) - Ephi(ir, iz)) / dz
				drEphi_dr = (Ephi(ir+1, iz) + Ephi(ir, iz))/(2.0d0*(real(ir, 8)+0.5d0)*dr) + dEphi_dr

				Ez(ir, iz) = (1-2*pi*dt*sigma_0(ir, iz))* Ez(ir, iz) /(1+2*pi*dt*sigma_0(ir, iz))  -&
							 4.0d0 * pi * dt * jz(ir, iz) /(1+2*pi*dt*sigma_0(ir, iz)) +&
							  c * dt * drBphi_dr/(1+2*pi*dt*sigma_0(ir, iz)) 
				Bz(ir, iz) = Bz(ir, iz) - c * dt * drEphi_dr
				end do				
			end do

	!$omp end parallel do 

	!$omp parallel do private(dBphi_dz, dEphi_dz, iz)
        do iz=1, nz-2
			do ir = 1, nr-2
				dBphi_dz= (Bphi(ir, iz+1) - Bphi(ir,  iz)) /  dz
				
				dEphi_dz = (Ephi(ir, iz+1) - Ephi(ir, iz)) / dz

				Er(ir, iz) = (1-2*pi*dt*sigma_0(ir, iz))* Er(ir, iz) /(1+2*pi*dt*sigma_0(ir, iz)) - &
							4.0d0 * pi * dt *jr(ir, iz)/(1+2*pi*dt*sigma_0(ir, iz))  - &
							c * dt * dBphi_dz/(1+2*pi*dt*sigma_0(ir, iz)) 
	
				Br(ir, iz) = Br(ir, iz) + c * dt * dEphi_dz
	
				jphi(ir, iz) = jphi(ir, iz) * (1 - dt * nu(ir,iz)/2)/ (1 + dt * nu(ir,iz)/2)+ &
								dt * omega_p2(ir, iz) / (4.0d0 * pi) * Ephi(ir, iz)/ (1 + dt * nu(ir,iz)/2) + &
									dt * omega_ce  *  jr(ir, iz)/ (1 + dt * nu(ir,iz)/2)
				end do				
			end do

	!$omp end parallel do 
	
	
    !$omp parallel do private(iz, dBr_dz, dBz_dr, dEr_dz, dEz_dr)
		do iz=1, nz-2
			do ir = 1, nr-2
				dBr_dz = (Br(ir, iz)- Br(ir, iz-1) ) / dz
				dBz_dr = (Bz(ir, iz) - Bz(ir-1, iz) ) / dr

				dEr_dz = ( Er(ir, iz) - Er(ir, iz-1)) / dz
				dEz_dr = (Ez(ir, iz)  - Ez(ir-1, iz)) / dr

				jr(ir, iz) = jr(ir, iz) * (1 - dt * nu(ir,iz))/ (1 + dt * nu(ir,iz)/2) + &
							dt * omega_p2(ir, iz)/ (4.0d0 * pi) *  Er(ir, iz)/ (1 + dt * nu(ir,iz)/2) - &
							dt * omega_ce *  jphi(ir, iz)/ (1 + dt * nu(ir,iz)/2)
	
				Ephi(ir, iz) =  (1-2*pi*dt*sigma_0(ir, iz))* Ephi(ir, iz) /(1+2*pi*dt*sigma_0(ir, iz)) - &
							4.0d0 * pi * dt * jphi(ir, iz) /(1+2*pi*dt*sigma_0(ir, iz))+ &
							c * dt * (dBr_dz-dBz_dr)/(1+2*pi*dt*sigma_0(ir, iz))
		
				Bphi(ir, iz) = Bphi(ir, iz) - c * dt * dEr_dz + c * dt * dEz_dr
			end do
		end do
	!$omp end parallel do   	
	
	!$omp parallel do
		do iz=0, nz-2
			do ir = 1, nr-2	
				jz(ir, iz) = jz(ir, iz)* (1 - dt * nu(ir, iz))/ (1 + dt * nu(ir,iz)/2) +&
							 dt * omega_p2(ir,iz)/ (4.0d0 * pi) *  Ez(ir, iz)/ (1 + dt * nu(ir,iz)/2)
			end do
		end do
	!$omp end parallel do 

		View(it) = Ephi(1, nint(real(nz, 8)/2))
	end do

	!call cpu_time(stopTime)
    call system_clock(int2)
 	 print *, 'Elapsed time, m : ',  (int2-int1)/1000/60 !((stopTime - startTime)/60.0d0)

	!open(24, file = 'View.bin', action = 'write', form = 'unformatted')
	!	write(24) nt
	!	write(24) dt
	!	write(24) View
	!close(24)

	!print * , max(real(jz), imag(jz))
	open(22, file = 'duct.bin', action = 'write', form = 'unformatted')
		write(22) nr
		write(22) nz
		write(22) dr
		write(22) dz
		!write(22) jext
		write(22) Er
		write(22) Ephi
		write(22) Ez
		write(22) Br
		write(22) Bphi
		write(22) Bz
	close(22)

	
	
	

	deallocate(Er, Ephi, Ez, Br, Bphi, Bz, jext)
	deallocate(Ne, nu)
	deallocate(jr, jphi, jz, sigma_0)
	deallocate(omega_p2, View)

	 
contains

	subroutine source(f, dt, dr, dz, Ra, t, r_source, z_source, jext)
		real(8), intent (in) :: f, dt, dr, dz, Ra
		integer(4), intent (in) :: t, r_source, z_source
		real(8) :: tau
		real(8) :: R1, R2, sigmaR, sigmaZ, k0
		real(8), dimension(:, :), intent(inout) :: jext

		jext(:,:) = (0d0, 0d0)
		tau = 5.0d0 / f 
		sigmaR = 0.5d0
		sigmaZ = 0.2d0 !разделить экспоненту R1**2/(2*sigmaR**2) с разными сигма
		k0 = 2*pi*f/c
		do iz = 1, nz
			do ir = 1, nr
			R1 = ((real(ir-r_source, 8)*dr)**2+(real(iz-z_source, 8)*dz)**2)**0.5
				if (R1<sigmaR) then
					jext(ir, iz) =  tanh(t * dt / tau) * sin(2.0d0 * pi * f * t* dt) * &
					exp(-R1**2/(2*sigmaR**2))!*exp(-j*k0 * iz*dz)
				end if
			end do
		end do
	end subroutine

	subroutine concentration(Ne0, Ne)
		real(8), dimension(:,:), intent (in) :: Ne0
		integer(4) :: i, j
		real(8) :: R, k, dn1, dn2, z1, z2
		real(8), dimension(:,:), intent (out) :: Ne
		
		R = 20.0d0
		dn1 = -0.2d0
		dn2 = 0.02d0
		z1 = 45.0d0
		z2 = nz*dz-40.0d0!z1
		k = (dn2-dn1)/(z2 - z1)
		!$omp parallel do
	    do j = LBOUND(Ne, 2), UBOUND(Ne, 2)
       		do i = LBOUND(Ne, 1), UBOUND(Ne, 1)
				if (i*dr<R) then
					Ne(i, j) = Ne0(i, j) *(1+dn1*((1+cos(pi*i*dr/(R)))/2.0d0))
					!if (j*dz<z1)  then
					!	Ne(i, j) = Ne0(i, j) *(1+dn1*((1+cos(pi*i*dr/(R)))/2.0d0))
					!else if (j*dz>z2) then
					!	Ne(i, j) = Ne0(i, j) *(1+dn2*((1+cos(pi*i*dr/(R)))/2.0d0))
					!else
					!	Ne(i, j) = Ne0(i, j) *(1+(k * (j * dz - z1)+ dn1)*((1+cos(pi*i*dr/(R)))/2.0d0))
					!end if
				end if
        	end do
		end do
		!$omp end parallel do
	end subroutine 

	subroutine absorb_layer(f, sigma)
		real(8), intent (in) :: f
		integer(4) :: isigma, i, j, jsigma
		real(8) :: sigmaMAX, lambda0
		real(8), dimension(:,:), intent (inout) :: sigma
		lambda0 = 10.0d0!c/f
		isigma = 3*nint(lambda0/dr)
		jsigma = 6*nint(lambda0/dz)
		sigmaMAX = 1.0d12/9.0d9
		
		!$omp parallel do
		do j = LBOUND(sigma, 2), UBOUND(sigma, 2)
	   		do i = LBOUND(sigma, 1), isigma
       	        sigma(nr-i+1, j) = sigmaMAX*((isigma-i)**3) ! 1d-6*
        	end do
		end do
		!$omp end parallel do

		!$omp parallel do
		do j = LBOUND(sigma, 2), jsigma
			do i = LBOUND(sigma, 1), UBOUND(sigma, 1)
    	        !sigma(i, j) =  max(sigmaMAX*(jsigma-j)**3, sigma(i, j))!1d-6*
    	    	sigma(i, nz-j+1) = max(sigmaMAX*(jsigma-j)**3, sigma(i, j))!sigma(i, j)
    	   end do
    	end do
		!$omp end parallel do
	end subroutine 

	
    !omega_p(Ne)**2
    subroutine omega_psq(Ne, omega_p2)
        real(8), dimension(:,:), intent(in) :: Ne
		real(8), dimension(:,:), intent(out) :: omega_p2
        omega_p2 = Ne * 4.0d0 * pi * e**2 / me
    end subroutine 
    
    !omega_ce(B)
    subroutine omega_c(B, omega_ce)
        real(8), intent(in) :: B
		real(8), intent(out) :: omega_ce
		omega_ce = e * B / me / c
    end subroutine 
end program