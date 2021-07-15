! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE deck

  ! Basic operations
  USE strings
  ! Deck internals
  USE deck_constant_block
  ! Deck Blocks
  USE deck_control_block
  USE deck_stencil_block
  USE deck_boundaries_block
  USE deck_species_block
  USE deck_injector_block
  USE deck_io_block
  USE deck_io_global_block
  USE deck_window_block
  USE deck_subset_block
  USE deck_collision_block
  USE deck_part_from_file_block
#ifdef PHOTONS
  USE photons
#endif
#ifdef BREMSSTRAHLUNG
  USE bremsstrahlung
#endif
  USE deck_qed_block
  USE deck_bremsstrahlung_block
  ! Initial Condition Blocks
  USE deck_laser_block
  USE deck_fields_block
  ! Extended IO Blocks
  USE deck_dist_fn_block
#ifndef NO_PARTICLE_PROBES
  USE deck_particle_probe_block
#endif
  ! Custom blocks
  USE custom_deck
  USE utilities
  USE version_data
  USE sdf

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: read_deck, write_input_decks, deallocate_input_deck_buffer

  SAVE
  CHARACTER(LEN=string_length) :: current_block_name
  LOGICAL :: invalid_block

  INTEGER, PARAMETER :: buffer_size = 1024
  TYPE :: file_buffer
    CHARACTER(LEN=c_max_path_length) :: filename
    CHARACTER(LEN=buffer_size), DIMENSION(:), POINTER :: buffer
    CHARACTER(LEN=32) :: md5sum
    INTEGER :: pos, idx, length
    TYPE(file_buffer), POINTER :: next
  END TYPE file_buffer

  TYPE(file_buffer), POINTER :: file_buffer_head
  INTEGER :: nbuffers = 0
  LOGICAL :: has_errors

CONTAINS

  !----------------------------------------------------------------------------
  ! These subroutines actually call the routines which read the deck blocks
  !----------------------------------------------------------------------------

  SUBROUTINE deck_initialise

    has_errors = .FALSE.

    CALL boundary_deck_initialise
    CALL collision_deck_initialise
    CALL constant_deck_initialise
    CALL control_deck_initialise
    CALL dist_fn_deck_initialise
    CALL fields_deck_initialise
    CALL injector_deck_initialise
    CALL io_deck_initialise
    CALL io_global_deck_initialise
    CALL laser_deck_initialise
    CALL stencil_deck_initialise
    CALL subset_deck_initialise
#ifndef NO_PARTICLE_PROBES
    CALL probe_deck_initialise
#endif
    CALL qed_deck_initialise
    CALL bremsstrahlung_deck_initialise
    CALL species_deck_initialise
    CALL window_deck_initialise
    CALL part_from_file_deck_initialise

  END SUBROUTINE deck_initialise



  SUBROUTINE deck_finalise(errcode_deck)

    INTEGER, INTENT(INOUT) :: errcode_deck

    CALL check_compulsory_blocks(errcode_deck)
    CALL boundary_deck_finalise
    CALL collision_deck_finalise
    CALL constant_deck_finalise
    CALL control_deck_finalise
    CALL dist_fn_deck_finalise
    CALL fields_deck_finalise
    CALL injector_deck_finalise
    CALL io_deck_finalise
    CALL io_global_deck_finalise
    CALL laser_deck_finalise
    CALL stencil_deck_finalise
    CALL subset_deck_finalise
#ifndef NO_PARTICLE_PROBES
    CALL probe_deck_finalise
#endif
    CALL qed_deck_finalise
    CALL bremsstrahlung_deck_finalise
    CALL species_deck_finalise
    CALL part_from_file_deck_finalise ! Must be called after
                                      ! species_deck_finalise
    CALL window_deck_finalise

  END SUBROUTINE deck_finalise


  ! This subroutine is called when a new block is started
  ! If a block NEEDS to do something when it starts, then
  ! The revelevant subroutine should be called here
  SUBROUTINE start_block(block_name)

    CHARACTER(LEN=*), INTENT(IN) :: block_name

    IF (str_cmp(block_name, 'boundaries')) THEN
      CALL boundary_block_start
    ELSE IF (str_cmp(block_name, 'collisions')) THEN
      CALL collision_block_start
    ELSE IF (str_cmp(block_name, 'constant')) THEN
      CALL constant_block_start
    ELSE IF (str_cmp(block_name, 'control')) THEN
      CALL control_block_start
    ELSE IF (str_cmp(block_name, 'dist_fn')) THEN
      CALL dist_fn_block_start
    ELSE IF (str_cmp(block_name, 'fields')) THEN
      CALL fields_block_start
    ELSE IF (str_cmp(block_name, 'output')) THEN
      CALL io_block_start
    ELSE IF (str_cmp(block_name, 'output_global')) THEN
      CALL io_global_block_start
    ELSE IF (str_cmp(block_name, 'laser')) THEN
      CALL laser_block_start
    ELSE IF (str_cmp(block_name, 'injector')) THEN
      CALL injector_block_start
    ELSE IF (str_cmp(block_name, 'stencil')) THEN
      CALL stencil_block_start
    ELSE IF (str_cmp(block_name, 'subset')) THEN
      CALL subset_block_start
#ifndef NO_PARTICLE_PROBES
    ELSE IF (str_cmp(block_name, 'probe')) THEN
      CALL probe_block_start
#endif
    ELSE IF (str_cmp(block_name, 'qed')) THEN
      CALL qed_block_start
    ELSE IF (str_cmp(block_name, 'bremsstrahlung')) THEN
      CALL bremsstrahlung_block_start
    ELSE IF (str_cmp(block_name, 'species')) THEN
      CALL species_block_start
    ELSE IF (str_cmp(block_name, 'window')) THEN
      CALL window_block_start
    ELSE IF (str_cmp(block_name, 'particles_from_file')) THEN
      CALL part_from_file_block_start
    END IF

  END SUBROUTINE start_block



  ! This subroutine is called when a new block is ended
  ! If a block NEEDS to do something when it ends, then
  ! The revelevant subroutine should be called here
  SUBROUTINE end_block(block_name)

    CHARACTER(LEN=*), INTENT(IN) :: block_name

    IF (str_cmp(block_name, 'boundaries')) THEN
      CALL boundary_block_end
    ELSE IF (str_cmp(block_name, 'collisions')) THEN
      CALL collision_block_end
    ELSE IF (str_cmp(block_name, 'constant')) THEN
      CALL constant_block_end
    ELSE IF (str_cmp(block_name, 'control')) THEN
      CALL control_block_end
    ELSE IF (str_cmp(block_name, 'dist_fn')) THEN
      CALL dist_fn_block_end
    ELSE IF (str_cmp(block_name, 'fields')) THEN
      CALL fields_block_end
    ELSE IF (str_cmp(block_name, 'output')) THEN
      CALL io_block_end
    ELSE IF (str_cmp(block_name, 'output_global')) THEN
      CALL io_global_block_end
    ELSE IF (str_cmp(block_name, 'laser')) THEN
      CALL laser_block_end
    ELSE IF (str_cmp(block_name, 'injector')) THEN
      CALL injector_block_end
    ELSE IF (str_cmp(block_name, 'stencil')) THEN
      CALL stencil_block_end
    ELSE IF (str_cmp(block_name, 'subset')) THEN
      CALL subset_block_end
#ifndef NO_PARTICLE_PROBES
    ELSE IF (str_cmp(block_name, 'probe')) THEN
      CALL probe_block_end
#endif
    ELSE IF (str_cmp(block_name, 'qed')) THEN
      CALL qed_block_end
    ELSE IF (str_cmp(block_name, 'bremsstrahlung')) THEN
      CALL bremsstrahlung_block_end
    ELSE IF (str_cmp(block_name, 'species')) THEN
      CALL species_block_end
    ELSE IF (str_cmp(block_name, 'window')) THEN
      CALL window_block_end
    ELSE IF (str_cmp(block_name, 'particles_from_file')) THEN
      CALL part_from_file_block_end
    END IF

  END SUBROUTINE end_block



  FUNCTION handle_block(block_name, block_element, block_value)

    CHARACTER(LEN=*), INTENT(IN) :: block_name, block_element, block_value
    INTEGER :: handle_block, io, iu
    LOGICAL, SAVE :: deo_warn = .TRUE.

    handle_block = c_err_none
    ! Constants can be defined in any deck state, so put them here
    IF (str_cmp(block_name, 'constant') &
            .OR. str_cmp(block_name, 'deo')) THEN
      IF (rank == 0 .AND. str_cmp(block_name, 'deo') .AND. deo_warn) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'The block name "deo" is deprecated.'
          WRITE(io,*) 'Please use the block name "constant" instead.'
        END DO
        deo_warn = .FALSE.
      END IF
      handle_block = constant_block_handle_element(block_element, block_value)
      RETURN
    END IF

    ! Test for known blocks
    IF (str_cmp(block_name, 'boundaries')) THEN
      handle_block = boundary_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'collisions')) THEN
      handle_block = collision_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'constant')) THEN
      handle_block = constant_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'control')) THEN
      handle_block = control_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'dist_fn')) THEN
      handle_block = dist_fn_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'fields')) THEN
      handle_block = fields_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'output')) THEN
      handle_block = io_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'output_global')) THEN
      handle_block = io_global_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'laser')) THEN
      handle_block = laser_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'injector')) THEN
      handle_block = injector_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'stencil')) THEN
      handle_block = stencil_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'subset')) THEN
      handle_block = subset_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'probe')) THEN
#ifndef NO_PARTICLE_PROBES
      handle_block = probe_block_handle_element(block_element, block_value)
#else
      handle_block = IOR(handle_block, c_err_pp_options_wrong)
      extended_error_string = '-DNO_PARTICLE_PROBES'
#endif
      RETURN
    ELSE IF (str_cmp(block_name, 'qed')) THEN
      handle_block = qed_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'bremsstrahlung')) THEN
      handle_block = bremsstrahlung_block_handle_element(block_element, &
          block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'species')) THEN
      handle_block = species_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'window')) THEN
      handle_block = window_block_handle_element(block_element, block_value)
      RETURN
    ELSE IF (str_cmp(block_name, 'particles_from_file')) THEN
      handle_block = &
          part_from_file_block_handle_element(block_element, block_value)
      RETURN
    END IF

    handle_block = c_err_unknown_block
    ! Pass through to the custom block
    handle_block = custom_blocks_handle_element(block_name, block_element, &
        block_value)

  END FUNCTION handle_block



  ! These subroutines are there to check for the basic minimal compulsory
  ! blocks are present. They're a bit ugly, but they seem to be the easiest
  ! way to do it without adding complexity to the code
  SUBROUTINE check_compulsory_blocks(errcode_deck)

    LOGICAL :: problem_found
    INTEGER, INTENT(INOUT) :: errcode_deck
    INTEGER :: io, iu

    IF (deck_state == c_ds_first) RETURN

    problem_found = .FALSE.

    errcode_deck = c_err_none
    errcode_deck = IOR(errcode_deck, boundary_block_check())
    IF (use_qed) THEN
#ifdef PHOTONS
      errcode_deck = IOR(errcode_deck, check_qed_variables())
#endif
      errcode_deck = IOR(errcode_deck, qed_block_check())
    END IF
    IF (use_bremsstrahlung) THEN
#ifdef BREMSSTRAHLUNG
      errcode_deck = IOR(errcode_deck, check_bremsstrahlung_variables())
#endif
      errcode_deck = IOR(errcode_deck, bremsstrahlung_block_check())
    END IF
    errcode_deck = IOR(errcode_deck, constant_block_check())
    errcode_deck = IOR(errcode_deck, control_block_check())
    errcode_deck = IOR(errcode_deck, dist_fn_block_check())
    errcode_deck = IOR(errcode_deck, fields_block_check())
    errcode_deck = IOR(errcode_deck, io_block_check())
    errcode_deck = IOR(errcode_deck, io_global_block_check())
    errcode_deck = IOR(errcode_deck, laser_block_check())
    errcode_deck = IOR(errcode_deck, injector_block_check())
    errcode_deck = IOR(errcode_deck, stencil_block_check())
    errcode_deck = IOR(errcode_deck, subset_block_check())
#ifndef NO_PARTICLE_PROBES
    errcode_deck = IOR(errcode_deck, probe_block_check())
#endif
    errcode_deck = IOR(errcode_deck, species_block_check())
    errcode_deck = IOR(errcode_deck, window_block_check())
    errcode_deck = IOR(errcode_deck, part_from_file_block_check())

    errcode_deck = IOR(errcode_deck, custom_blocks_check())

    problem_found = (IAND(errcode_deck, c_err_missing_elements) /= 0)

    IF (problem_found) THEN
      errcode_deck = IOR(errcode_deck, c_err_terminate)
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Not all required elements of input deck specified.'
          WRITE(io,*) 'Please fix input deck and rerun code'
        END DO
      END IF
    ELSE
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) 'Initial conditions complete and valid.'
        END DO
      END IF
    END IF

  END SUBROUTINE check_compulsory_blocks



  !----------------------------------------------------------------------------
  ! These subroutines are the in depth detail of how the parser works
  !----------------------------------------------------------------------------

  RECURSIVE SUBROUTINE read_deck(filename, first_call, deck_state_in)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    LOGICAL, INTENT(IN) :: first_call
    INTEGER, INTENT(IN) :: deck_state_in
    CHARACTER :: u0, u1
    INTEGER :: pos = 1, flip = 1, slen, s, f, elements = 0, lun
    INTEGER :: line, column
    LOGICAL :: ignore, continuation
    LOGICAL, SAVE :: warn = .TRUE.
    LOGICAL, SAVE :: warn_ascii = .FALSE.
    TYPE(string_type), DIMENSION(2) :: deck_values
    CHARACTER(LEN=c_max_path_length) :: deck_filename
    CHARACTER(LEN=c_max_path_length) :: const_filename
    CHARACTER(LEN=string_length) :: len_string
    LOGICAL :: terminate = .FALSE.
    LOGICAL :: exists
    INTEGER :: errcode_deck, i, io, iu, rank_check
    CHARACTER(LEN=buffer_size), DIMENSION(:), ALLOCATABLE :: tmp_buffer
    TYPE(file_buffer), POINTER :: fbuf
    LOGICAL :: already_parsed, got_eor, got_eof

    deck_state = deck_state_in
    ! No error yet
    errcode_deck = c_err_none
    ! Characteristic string which represents a 'blank' string
    blank = 'BLANKBLANK'

    lun = 5
    rank_check = 0
    already_parsed = .FALSE.
    u0 = ' '

    ! Make the whole filename by adding the data_dir to the filename
    deck_filename = TRIM(ADJUSTL(data_dir)) // '/' // TRIM(ADJUSTL(filename))

    ! deck_state tells the code whether it's parsing the normal input deck
    ! Or the initial conditions. You can add more states if you want.
    ! Just search for deck_state

    ! If this is the first time that this deck has been called then do some
    ! housekeeping. Put any initialisation code that is needed in here
    IF (first_call) CALL deck_initialise

#ifdef DECK_DEBUG
    IF (deck_state == c_ds_first) all_deck_errcodes_fatal = .TRUE.
#endif

    ! Flag which tells the code when a # or \ character has been
    ! found and everything beyond it is to be ignored
    ignore = .FALSE.
    continuation = .FALSE.
    status_filename = TRIM(ADJUSTL(data_dir)) // '/deck.status'
    const_filename = TRIM(ADJUSTL(data_dir)) // '/const.status'

    ! rank 0 reads the file and then passes it out to the other nodes using
    ! MPI_BCAST
    IF (rank == 0) THEN
      IF (.NOT. ASSOCIATED(file_buffer_head)) THEN
        ALLOCATE(file_buffer_head)
        fbuf => file_buffer_head
        fbuf%filename = ''
        NULLIFY(fbuf%next)
        NULLIFY(fbuf%buffer)
      ELSE
        fbuf => file_buffer_head
      END IF

      DO WHILE (ASSOCIATED(fbuf%next))
        fbuf => fbuf%next
        IF (fbuf%filename == deck_filename) THEN
          already_parsed = .TRUE.
          EXIT
        END IF
      END DO
      IF (.NOT. already_parsed) THEN
        ALLOCATE(fbuf%next)
        nbuffers = nbuffers + 1
        fbuf => fbuf%next
        fbuf%filename = deck_filename
        fbuf%pos = 1
        fbuf%idx = 1
        fbuf%length = 1
        NULLIFY(fbuf%next)
        ALLOCATE(fbuf%buffer(fbuf%length))
      END IF

      ! Check whether or not the input deck file requested exists
      INQUIRE(file=deck_filename, exist=exists)
      IF (.NOT. exists .AND. rank == 0) THEN
        PRINT *, '*** ERROR ***'
        PRINT *, 'Input deck file "' // TRIM(deck_filename) &
            // '" does not exist.'
        PRINT *, 'Create the file and rerun the code.'
        CALL abort_code(c_err_io_error)
      END IF

      ! Get a free lun. Don't use a constant lun to allow for recursion
      lun = get_free_lun()
      OPEN(unit=lun, file=TRIM(ADJUSTL(deck_filename)))
#ifndef NO_IO
      IF (first_call .AND. rank == 0) THEN
        ! Create a new file on first pass, otherwise append
        IF (deck_state == c_ds_first) THEN
          OPEN(unit=du, status='REPLACE', file=status_filename, iostat=errcode)
          WRITE(du,*) ascii_header
          WRITE(du,*)
        ELSE
          OPEN(unit=du, status='OLD', position='APPEND', file=status_filename, &
              iostat=errcode)
        END IF
        IF (print_deck_constants) THEN
          OPEN(unit=duc, status='REPLACE', file=const_filename, iostat=errcode)
        END IF

        WRITE(du,'(a,i3)') 'Deck state:', deck_state
        WRITE(du,*)
      END IF
#endif
      deck_values(1)%value = ''
      deck_values(2)%value = ''
      slen = 1
      line = 1
      deck_line_number = '0'
      column = 0

      ! Use non-advancing IO to pop characters off the deck file one at a time
      ! Use basic token parsing to split into two substrings across
      ! an '=' or ':' symbol
      DO
        errcode_deck = c_err_none
        ! Read a character
        ! For ordinary characters, f is zero
        ! When an end of line character is read then got_eor is .TRUE.
        ! When end of file is reached, f is negative and got_eor is .FALSE.
        got_eor = .TRUE.
        got_eof = .FALSE.
        READ(lun, '(A1)', advance='no', size=s, iostat=f, eor=10) u1

        IF (f < 0) THEN
          got_eor = .TRUE.
          got_eof = .TRUE.
          column = 0
          line = line + 1
          CALL integer_as_string(line - 1, deck_line_number)
        ELSE
          got_eor = .FALSE.
          column = column + 1
        END IF

        IF (.NOT.ignore .AND. IACHAR(u1) > 127) THEN
          IF (rank == rank_check) THEN
            WRITE(du, '(A,I5,A,I5,A,I5)') 'Line', line, ', column', column, &
                ', non-ascii character discarded:', IACHAR(u1)
          END IF
          warn_ascii = .TRUE.
          CYCLE
        END IF

10      IF (.NOT. already_parsed) THEN
          ! Store character in a buffer so that we can write the input deck
          ! contents to a restart dump
          IF (f == 0) THEN
            fbuf%buffer(fbuf%idx)(fbuf%pos:fbuf%pos) = u1
          ELSE IF (got_eor .AND. .NOT.got_eof) THEN
            fbuf%buffer(fbuf%idx)(fbuf%pos:fbuf%pos) = ACHAR(10) ! new line
          ELSE
            fbuf%buffer(fbuf%idx)(fbuf%pos:fbuf%pos) = ACHAR(0)  ! null
            fbuf%pos = fbuf%pos - 1
          END IF

          ! If we reached the end of the character string then move to the next
          ! element of the array
          IF (fbuf%pos == buffer_size) THEN
            ! If we reached the end of the array then allocate some more
            IF (fbuf%idx == fbuf%length) THEN
              ALLOCATE(tmp_buffer(fbuf%length))
              DO i = 1,fbuf%length
                tmp_buffer(i) = fbuf%buffer(i)
              END DO
              DEALLOCATE(fbuf%buffer)
              ALLOCATE(fbuf%buffer(2*fbuf%length))
              DO i = 1,fbuf%length
                fbuf%buffer(i) = tmp_buffer(i)
              END DO
              DEALLOCATE(tmp_buffer)
              fbuf%length = 2*fbuf%length
            END IF
            fbuf%pos = 1
            fbuf%idx = fbuf%idx + 1
          ELSE
            fbuf%pos = fbuf%pos + 1
          END IF
        END IF

        IF (continuation .AND. warn) THEN
          IF (u1 /= ' ' .AND. u1 /= ACHAR(9)) THEN ! ACHAR(9) = tab
            IF (rank == rank_check) THEN
              DO iu = 1, nio_units ! Print to stdout and to file
                io = io_units(iu)
                WRITE(io,*)
                WRITE(io,*) '*** WARNING ***'
                WRITE(io,*) 'Extra characters after continuation line in', &
                    ' input deck.'
                WRITE(io,*)
              END DO
            END IF
            warn = .FALSE.
          END IF
        END IF

        ! If the character is a # or \ then ignore the rest of the line
        IF (u1 == '#') THEN
          ignore = .TRUE.
        ELSE IF (u1 == ACHAR(92)) THEN ! ACHAR(92) = '\'
          ignore = .TRUE.
          continuation = .TRUE.
        END IF

        ! If not in comment mode then use the character
        IF (.NOT. ignore) THEN
          ! If the current character isn't a special character then just stick
          ! it in the buffer
          ! ACHAR(9) = tab
          IF (u1 /= '=' .AND. u1 /= ACHAR(9) .AND. u1 /= ':' &
              .AND. u1 /= ACHAR(13) .AND. f == 0) THEN
            IF (u1 /= ' ' .OR. u0 /= ' ') THEN
              deck_values(flip)%value(pos:pos) = u1
              pos = pos + 1
              slen = slen + 1
              u0 = u1
              IF (pos > string_length) pos = string_length
            END IF
          END IF

          ! If it's equals or : then you're parsing the other part of the
          ! expression
          IF (u1 == '=' .OR. u1 == ':') THEN
            flip = 2
            pos = 1
            slen = 1
          END IF
        END IF

        ! If got_eor is .TRUE. then you've reached the end of the line, so
        ! reset comment and continuation states
        IF (got_eor) THEN
          ignore = .FALSE.
          IF (continuation) THEN
            got_eor = .FALSE.
            f = 0
          END IF
          continuation = .FALSE.
          line = line + 1
          CALL integer_as_string(line - 1, deck_line_number)
          column = 0
        END IF

        ! If you've not read a blank line then
        IF (got_eor .AND. pos > 1) THEN
          IF (slen > string_length) THEN
            CALL integer_as_string(slen, len_string)
            DO iu = 1, nio_units ! Print to stdout and to file
              io = io_units(iu)
              WRITE(io,*)
              WRITE(io,*) '*** ERROR ***'
              WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
              IF (flip > 1) THEN
                WRITE(io,*) 'Whilst reading ',TRIM(deck_values(1)%value) &
                    // ' = ' // TRIM(deck_values(2)%value(1:pos-1))
              ELSE
                WRITE(io,*) 'Whilst reading ', &
                    TRIM(deck_values(1)%value(1:pos-1))
              END IF
              WRITE(io,*) 'String value too long. Please increase the size ', &
                  'of "string_length" in ','shared_data.F90 to be at least ', &
                  TRIM(len_string)
            END DO
            CALL abort_code(c_err_io_error)
          END IF
          elements = elements+1
          flip = 1
          pos = 1
          slen = 1
          deck_values(1)%value = TRIM(ADJUSTL(deck_values(1)%value))
          deck_values(2)%value = TRIM(ADJUSTL(deck_values(2)%value))
          CALL MPI_BCAST(1, 1, MPI_INTEGER, 0, comm, errcode)
          CALL MPI_BCAST(nbuffers, 1, MPI_INTEGER, 0, comm, errcode)
          CALL MPI_BCAST(deck_values(1)%value, string_length, MPI_CHARACTER, &
              0, comm, errcode)
          CALL MPI_BCAST(deck_values(2)%value, string_length, MPI_CHARACTER, &
              0, comm, errcode)
          CALL handle_deck_element(deck_values(1)%value, deck_values(2)%value, &
              errcode_deck)
          deck_values(1)%value = ''
          deck_values(2)%value = ''
          ignore = .FALSE.
          continuation = .FALSE.
          u0 = ' '
        END IF
        IF (errcode_deck /= c_err_none) has_errors = .TRUE.
        IF (got_eof) THEN
          CALL MPI_BCAST(0, 1, MPI_INTEGER, 0, comm, errcode)
          CLOSE(lun)
          EXIT
        END IF
        terminate = terminate .OR. IAND(errcode_deck, c_err_terminate) /= 0
        IF (terminate) EXIT
      END DO

      IF (.NOT. already_parsed) THEN
        CALL sdf_md5_init()
        DO i = 1, fbuf%idx - 1
           fbuf%md5sum = sdf_md5_append(fbuf%buffer(i)(1:buffer_size))
        END DO
        fbuf%md5sum = sdf_md5_append(fbuf%buffer(fbuf%idx)(1:fbuf%pos-1))
        IF (MOD(fbuf%pos-1, 64) == 0) fbuf%md5sum = sdf_md5_append('')
      END IF
    ELSE
      DO
        errcode_deck = c_err_none
        CALL MPI_BCAST(f, 1, MPI_INTEGER, 0, comm, errcode)
        IF (f == 0) EXIT
          CALL MPI_BCAST(nbuffers, 1, MPI_INTEGER, 0, comm, errcode)
        CALL MPI_BCAST(deck_values(1)%value, string_length, MPI_CHARACTER, &
              0, comm, errcode)
        CALL MPI_BCAST(deck_values(2)%value, string_length, MPI_CHARACTER, &
              0, comm, errcode)
        CALL handle_deck_element(deck_values(1)%value, deck_values(2)%value, &
              errcode_deck)
        deck_values(1)%value = ''
        deck_values(2)%value = ''
        terminate = terminate .OR. IAND(errcode_deck, c_err_terminate) /= 0
        IF (terminate) EXIT
      END DO
    END IF

    CALL MPI_BARRIER(comm, errcode)
!!$    IF (.NOT. first_call)  THEN
!!$       RETURN
!!$    END IF

    ! Don't check compulsory blocks if going to bomb anyway, just stinks up
    ! the output file
    IF (.NOT. terminate .AND. first_call) CALL deck_finalise(errcode_deck)

    terminate = terminate .OR. IAND(errcode_deck, c_err_terminate) /= 0
    ! Fatal error, cause code to bomb
    IF (terminate .AND. rank == 0) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*)
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) 'The code cannot parse the input deck sufficiently to run.'
      END DO
      WRITE(*, *) 'Please read the output file "', TRIM(status_filename), &
          '" to check for errors.'
#ifndef NO_IO
      WRITE(du,*) 'Please read this file and correct any errors mentioned.'
#endif
    END IF

    IF (deck_state == c_ds_first .AND. warn_ascii .AND. rank == 0) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = io_units(iu)
        WRITE(io,*)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'Non-ascii characters present in the input deck ', &
            'were discarded'
      END DO
    END IF

#ifndef NO_IO
    IF (first_call .AND. rank == 0) THEN
      CLOSE(du)
      CLOSE(duc)
    END IF
#endif

    IF (terminate) CALL abort_code(c_err_generic_error)

    IF (all_deck_errcodes_fatal .AND. has_errors) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Deck has warnings and you have requested ' &
              // 'deck_warnings_fatal.'
          WRITE(io,*) 'Please fix input deck and rerun ' &
              // 'code or disable this option'
        END DO
      END IF
      CALL abort_code(c_err_generic_error)
    END IF

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE read_deck



  SUBROUTINE handle_deck_element(element, value, errcode_deck)

    CHARACTER(*), INTENT(IN) :: element
    CHARACTER(*), INTENT(IN) :: value
    INTEGER, INTENT(INOUT) :: errcode_deck
    INTEGER :: rank_check, io, iu
    INTEGER, SAVE :: err_count

    rank_check = 0

    IF (str_cmp(element, 'import')) THEN
#ifndef NO_IO
      IF (rank == rank_check) THEN
        WRITE(du,*)
        WRITE(du,*) 'Importing "' // TRIM(ADJUSTL(value)) // '" file'
        WRITE(du,*)
      END IF
#endif
      CALL read_deck(TRIM(ADJUSTL(value)), .FALSE., deck_state)
      RETURN
    END IF

    IF (str_cmp(element, 'begin')) THEN
      errcode_deck = handle_block(value, blank, blank)
      invalid_block = IAND(errcode_deck, c_err_unknown_block) /= 0
      invalid_block = invalid_block &
          .OR. (IAND(errcode_deck, c_err_pp_options_missing) /= 0) &
          .OR. (IAND(errcode_deck, c_err_pp_options_wrong) /= 0)
      IF (invalid_block .AND. rank == rank_check) THEN
        IF (IAND(errcode_deck, c_err_pp_options_missing) /= 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
            WRITE(io,*) 'The block "' // TRIM(value) &
                // '" cannot be set because'
            WRITE(io,*) 'the code has not been compiled with the correct ' &
                // 'preprocessor options.'
            WRITE(io,*) 'Code will continue, but to use selected features, ' &
                // 'please recompile with the'
            WRITE(io,*) TRIM(extended_error_string) // ' option'
            WRITE(io,*)
          END DO
        ELSE IF (IAND(errcode_deck, c_err_pp_options_wrong) /= 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
            WRITE(io,*) 'The block "' // TRIM(value) &
                // '" cannot be set because'
            WRITE(io,*) 'the code has not been compiled with the correct ' &
                // 'preprocessor options.'
            WRITE(io,*) 'Code will continue, but to use selected features, ' &
                // 'please recompile without the'
            WRITE(io,*) TRIM(extended_error_string) // ' option'
            WRITE(io,*)
          END DO
        ELSE
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
            WRITE(io,*) 'Unknown block "' // TRIM(value) &
                // '" in input deck, ignoring', deck_state
          END DO
        END IF
      END IF
      CALL start_block(value)
      err_count = 0
      current_block_name = value
#ifndef NO_IO
      IF (rank == rank_check) THEN
        WRITE(du,*) 'Beginning "' // TRIM(ADJUSTL(value)) // '" block'
        WRITE(du,*)
      END IF
#endif
      ! Reset errcode_deck here because reporting c_err_unknown_element is OK
      errcode_deck = c_err_none
      RETURN
    END IF
    IF (str_cmp(element, 'end')) THEN
      CALL end_block(current_block_name)
      invalid_block = .TRUE.
#ifndef NO_IO
      IF (rank == rank_check) THEN
        WRITE(du,*)
        WRITE(du,*) 'Ending "' // TRIM(ADJUSTL(value)) // '" block'
        WRITE(du,*)
        IF (err_count /= 0) THEN
          WRITE(du,*) '*** WARNING ***'
          WRITE(du,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(du,*) 'Block "' // TRIM(ADJUSTL(value)) // '" contains errors'
          WRITE(du,*)
        END IF
      END IF
#endif
      RETURN
    END IF

    ! Check invalid block to limit amount of rubbish that appears
    ! If the input deck is invalid
    IF (invalid_block) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Value "' // TRIM(value) // '" in element "' &
              // TRIM(element) // '" is invalid and cannot be parsed.'
          WRITE(io,*) 'Code will terminate'
          WRITE(io,*)
        END DO
      END IF
      errcode_deck = IOR(errcode_deck, c_err_terminate)
      RETURN
    ELSE
      errcode_deck = handle_block(current_block_name, element, value)
    END IF

#ifndef NO_IO
    IF (errcode_deck == c_err_none) THEN
      IF (rank == rank_check) &
          WRITE(du, *) ACHAR(9), 'Element ', TRIM(ADJUSTL(element)), '=', &
              TRIM(ADJUSTL(value)), ' handled OK'
      RETURN
    END IF
#endif
    ! Test for error conditions
    ! If an error is fatal then set terminate to .TRUE.
    IF (IAND(errcode_deck, c_err_unknown_element) /= 0) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Unrecognised element "' // TRIM(element) &
              // '" in input deck.'
          WRITE(io,*) 'Code will continue to run, but behaviour is undefined'
          WRITE(io,*)
        END DO
      END IF
    END IF
    IF (IAND(errcode_deck, c_err_preset_element) /= 0) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Element "' // TRIM(element) &
              // '" is set multiple times in this deck.'
          WRITE(io,*) 'Code will continue using first value in deck'
          WRITE(io,*)
        END DO
      END IF
    END IF
    IF (IAND(errcode_deck, c_err_preset_element_use_later) /= 0) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Element "' // TRIM(element) &
              // '" is set multiple times in this deck.'
          WRITE(io,*) 'Code will continue using last value in deck'
          WRITE(io,*)
        END DO
      END IF
    END IF
    IF (IAND(errcode_deck, c_err_bad_value) /= 0) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Value "' // TRIM(value) // '" in element "' &
              // TRIM(element) // '" is'
          WRITE(io,*) 'invalid or could not be parsed. Code will terminate.'
          WRITE(io,*)
        END DO
      END IF
      errcode_deck = IOR(errcode_deck, c_err_terminate)
    END IF
    IF (IAND(errcode_deck, c_err_warn_bad_value) /= 0) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Value "' // TRIM(value) // '" in element "' &
              // TRIM(element) // '" is'
          WRITE(io,*) 'invalid or could not be parsed. Code will use', &
              ' default value.'
          WRITE(io,*)
        END DO
      END IF
    END IF

    IF (IAND(errcode_deck, c_err_required_element_not_set) /= 0) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Value "' // TRIM(value) // '" in element "' &
              // TRIM(element) // '" cannot be'
          WRITE(io,*) 'set because a prerequisite element "' &
              // TRIM(extended_error_string) // '" has'
          WRITE(io,*) 'not been set. Code will terminate'
          WRITE(io,*)
        END DO
      END IF
      errcode_deck = IOR(errcode_deck, c_err_terminate)
    END IF
    IF (IAND(errcode_deck, c_err_pp_options_missing) /= 0) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'The element "' // TRIM(element) // '" of block "' &
              // TRIM(current_block_name) // '" cannot be set'
          WRITE(io,*) 'because the code has not been compiled with the ' &
              // 'correct preprocessor options.'
          WRITE(io,*) 'Code will continue, but to use selected features, ' &
              // 'please recompile with the'
          WRITE(io,*) TRIM(extended_error_string) // ' option'
          WRITE(io,*)
        END DO
      END IF
    END IF
    IF (IAND(errcode_deck, c_err_generic_warning) /= 0) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) TRIM(extended_error_string)
          WRITE(io,*)
        END DO
      END IF
    END IF
    IF (IAND(errcode_deck, c_err_generic_error) /= 0) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) TRIM(extended_error_string)
          WRITE(io,*)
        END DO
      END IF
      errcode_deck = IOR(errcode_deck, c_err_terminate)
    END IF
    IF (IAND(errcode_deck, c_err_other) /= 0) THEN
      IF (rank == rank_check) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'You have managed to find an impossible situation in ' &
              // 'this code.'
          WRITE(io,*) 'Good for you. Just because of that, the code will ' &
              // 'terminate.'
          WRITE(io,*)
        END DO
      END IF
      errcode_deck = IOR(errcode_deck, c_err_terminate)
    END IF

    err_count = err_count+1

  END SUBROUTINE handle_deck_element



  SUBROUTINE write_input_decks(handle)

    TYPE(sdf_file_handle) :: handle
    TYPE(file_buffer), POINTER :: fbuf
    CHARACTER(LEN=1) :: buffer(1)
    INTEGER :: i

    IF (rank == 0) THEN
      fbuf => file_buffer_head
      DO i = 1,nbuffers
        fbuf => fbuf%next

        CALL sdf_write_datablock(handle, 'input_deck/' // TRIM(fbuf%filename), &
            'EPOCH input deck: ' // TRIM(fbuf%filename), &
            fbuf%buffer(1:fbuf%idx-1), fbuf%buffer(fbuf%idx)(1:fbuf%pos-1), &
            'text/plain', 'md5', fbuf%md5sum)
      END DO
    ELSE
      ! These calls are required since sdf_write_datablock() is a collective
      ! operation
      DO i = 1,nbuffers
        CALL sdf_write_datablock(handle, 'input_deck/', 'EPOCH input deck: ', &
            buffer, buffer(1), 'text/plain', 'md5', '')
      END DO
    END IF

  END SUBROUTINE write_input_decks



  SUBROUTINE deallocate_input_deck_buffer

    TYPE(file_buffer), POINTER :: fbuf, next
    INTEGER :: stat

    IF (.NOT. ASSOCIATED(file_buffer_head)) RETURN

    fbuf => file_buffer_head%next
    DO WHILE (ASSOCIATED(fbuf))
      next => fbuf%next

      DEALLOCATE(fbuf%buffer, STAT=stat)
      DEALLOCATE(fbuf, STAT=stat)

      fbuf => next
    END DO

    DEALLOCATE(file_buffer_head, STAT=stat)
    nbuffers = 0

  END SUBROUTINE deallocate_input_deck_buffer

END MODULE deck
