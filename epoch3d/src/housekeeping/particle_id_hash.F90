MODULE particle_id_hash_mod

  USE constants
  USE shared_data
  USE random_generator
  USE utilities

  IMPLICIT NONE

  INTEGER, PARAMETER :: id_global_chunk_size = 1000000 ! 8MB chunks at i8
  LOGICAL :: random_state_set = .FALSE.
  TYPE(random_state_type), SAVE :: random_state

  TYPE :: particle_id_inner_list
    INTEGER(idkind), DIMENSION(:), ALLOCATABLE :: list
    CONTAINS
    PROCEDURE :: holds => pid_inner_list_holds
    PROCEDURE :: add => pid_inner_list_add
    PROCEDURE :: delete => pid_inner_list_delete
    FINAL :: pid_inner_list_destructor
  END TYPE particle_id_inner_list

  TYPE :: particle_id_hash
    PRIVATE
    CHARACTER(LEN=c_max_string_length) :: name
    TYPE(particle_id_inner_list), DIMENSION(:), ALLOCATABLE :: buckets
    INTEGER(i8) :: count, hash_gr
    INTEGER(i4) :: id_chunk_size
    CONTAINS
    PRIVATE
    PROCEDURE :: hash => pid_hash_hash
    PROCEDURE, PUBLIC :: holds => pid_hash_holds
    PROCEDURE, PUBLIC :: add => pid_hash_add
    PROCEDURE, PUBLIC :: add_if_local => pid_hash_add_if_local
    PROCEDURE, PUBLIC :: add_from_file => pid_hash_add_from_file
    PROCEDURE, PUBLIC :: delete => pid_hash_delete
    PROCEDURE :: init_i8 => pid_hash_init_i8
    PROCEDURE :: init_i4 => pid_hash_init_i4
    GENERIC, PUBLIC :: init => init_i8, init_i4
    PROCEDURE, PUBLIC :: optimise => pid_optimise
    FINAL :: pid_hash_destructor
  END TYPE particle_id_hash

  TYPE :: particle_id_hash_holder
    TYPE(particle_id_hash), POINTER :: contents => NULL()
  END TYPE particle_id_hash_holder

  TYPE :: particle_id_list_registry
    TYPE(particle_id_hash_holder), DIMENSION(:), ALLOCATABLE :: list
    CONTAINS
    PROCEDURE, PRIVATE :: add_new_hash => pidr_add_new_hash
    PROCEDURE, PRIVATE :: get_existing_hash_by_name => &
        pidr_get_existing_hash_by_name
    PROCEDURE, PRIVATE :: get_existing_hash_by_index => &
        pidr_get_existing_hash_by_index
    GENERIC :: get_existing_hash => get_existing_hash_by_name, &
        get_existing_hash_by_index
    PROCEDURE, PRIVATE :: get_hash_by_name => pidr_get_hash_by_name
    GENERIC :: get_hash => get_hash_by_name, get_existing_hash_by_index
    PROCEDURE :: get_hash_count => pidr_get_hash_count
    PROCEDURE :: map => pidr_map
    PROCEDURE :: delete_all => pidr_delete_all
    PROCEDURE :: delete_and_map => pidr_delete_and_map
    PROCEDURE :: add_with_map => pidr_add_with_map
    FINAL :: pidr_destructor
  END TYPE particle_id_list_registry

  TYPE(particle_id_list_registry), SAVE :: id_registry

  PRIVATE
  PUBLIC :: id_registry, particle_id_hash

CONTAINS

  !> Subroutine to sort an array in place

  RECURSIVE SUBROUTINE sort_list(list_in)

    INTEGER(idkind), DIMENSION(:), INTENT(INOUT) :: list_in
    INTEGER(i8) :: part_index

    IF (SIZE(list_in) <= 1) RETURN
    CALL partition(list_in, part_index)
    CALL sort_list(list_in(:part_index - 1))
    CALL sort_list(list_in(part_index:))

  END SUBROUTINE sort_list



  !> Function to partition an array on a pivot value
  !> Just use the first element as the pivot in the absence
  !> of a better scheme

  SUBROUTINE partition(list_in, part_index)

    INTEGER(idkind), DIMENSION(:), INTENT(INOUT) :: list_in
    INTEGER(idkind) :: pivot, temp, p1, p2, p3
    INTEGER(i8), INTENT(INOUT) :: part_index
    INTEGER(i8) :: upper, lower

    p1 = list_in(INT(random(random_state) * REAL(SIZE(list_in), num), idkind))
    p2 = list_in(INT(random(random_state) * REAL(SIZE(list_in), num), idkind))
    p3 = list_in(INT(random(random_state) * REAL(SIZE(list_in), num), idkind))

    pivot = median_of_three(p1, p2, p3)
    lower = 0
    upper = SIZE(list_in) + 1
    DO
      DO
        upper = upper - 1
        IF (list_in(upper) <= pivot) EXIT
      END DO
      DO
        lower = lower + 1
        IF (list_in(lower) >= pivot) EXIT
      END DO

      IF (lower < upper) THEN
        temp = list_in(lower)
        list_in(lower) = list_in(upper)
        list_in(upper) = temp
      ELSE IF (lower == upper) THEN
        part_index = lower + 1
        RETURN
      ELSE
        part_index = lower
        RETURN
      END IF
    END DO

  END SUBROUTINE partition



  FUNCTION median_of_three(a, b, c) RESULT(median)

    INTEGER(idkind), INTENT(IN) :: a, b, c
    INTEGER(idkind) :: median

    median = MAXVAL([MINVAL([a,b]), MINVAL([MAXVAL([a,b]), c])])

  END FUNCTION median_of_three



  !> Function to find if an item is in a list

  RECURSIVE FUNCTION in_list(list_in, test) RESULT(found)

    INTEGER(idkind), DIMENSION(:), INTENT(IN) :: list_in
    INTEGER(idkind), INTENT(IN) :: test
    LOGICAL :: found
    INTEGER(i8) :: trial_index

    IF (SIZE(list_in) == 0) THEN
      found = .FALSE.
      RETURN
    END IF

    IF (SIZE(list_in) == 1) THEN
      IF (list_in(1) == test) THEN
        found = .TRUE.
      ELSE
        found = .FALSE.
      END IF
      RETURN
    END IF

    trial_index = SIZE(list_in, KIND=i8) / 2_i8
    IF(test == list_in(trial_index)) THEN
      found = .TRUE.
      RETURN
    ELSE IF (test < list_in(trial_index)) THEN
      found = in_list(list_in(:trial_index-1), test)
    ELSE
      found = in_list(list_in(trial_index+1:), test)
    END IF

  END FUNCTION in_list



  !> Test if particle id is in this bin. Linear search

  FUNCTION pid_inner_list_holds(this, test_id, index_out) RESULT(holds)

    CLASS(particle_id_inner_list), INTENT(IN) :: this
    INTEGER(idkind), INTENT(IN) :: test_id
    INTEGER, INTENT(OUT), OPTIONAL :: index_out
    LOGICAL :: holds
    INTEGER :: ipart, sz, index_inner

    index_inner = -1
    holds = .FALSE.
    IF (ALLOCATED(this%list)) THEN
      sz = SIZE(this%list)
      DO ipart = 1, sz
        IF (this%list(ipart) == test_id) THEN
          holds = .TRUE.
          index_inner = ipart
          EXIT
        END IF
      END DO
    END IF
    IF (PRESENT(index_out)) index_out = index_inner

  END FUNCTION pid_inner_list_holds



  !> Add a particle to this list

  SUBROUTINE pid_inner_list_add(this, add_id)

    CLASS(particle_id_inner_list), INTENT(INOUT) :: this
    INTEGER(idkind), INTENT(IN) :: add_id
    INTEGER(idkind), DIMENSION(:), ALLOCATABLE :: temp
    LOGICAL :: copyback
    INTEGER :: sz

    copyback = ALLOCATED(this%list)
    sz = 0
    IF (copyback) THEN
      sz = SIZE(this%list)
      ALLOCATE(temp(sz), SOURCE = this%list)
      DEALLOCATE(this%list)
    END IF
    sz = sz + 1
    ALLOCATE(this%list(sz))
    this%list(sz) = add_id
    IF (copyback) THEN
      this%list(1:sz-1) = temp(1:sz-1)
      DEALLOCATE(temp)
    END IF

  END SUBROUTINE pid_inner_list_add



  !> Remove a particle from this list

  FUNCTION pid_inner_list_delete(this, del_id) RESULT(holds)

    CLASS(particle_id_inner_list), INTENT(INOUT) :: this
    INTEGER(idkind), INTENT(IN) :: del_id
    INTEGER(idkind), DIMENSION(:), ALLOCATABLE :: temp
    INTEGER :: sz, ind
    LOGICAL :: holds

    holds = this%holds(del_id, index_out = ind)

    IF (.NOT. holds) RETURN
    sz = SIZE(this%list)
    IF (sz == 1) THEN
      DEALLOCATE(this%list)
      RETURN
    END IF

    ALLOCATE(temp(sz-1))
    temp(1:ind - 1) = this%list(1:ind-1)
    temp(ind:sz-1) = this%list(ind+1:sz)
    DEALLOCATE(this%list)
    ALLOCATE(this%list(1:sz-1), SOURCE = temp)
    DEALLOCATE(temp)

  END FUNCTION pid_inner_list_delete



  !> Clean up inner list on destruction

  PURE ELEMENTAL SUBROUTINE pid_inner_list_destructor(this)

    TYPE(particle_id_inner_list), INTENT(INOUT) :: this

    IF (.NOT. ALLOCATED(this%list)) RETURN
    DEALLOCATE(this%list)

  END SUBROUTINE pid_inner_list_destructor



  !> Function to generate a hash from a particle id

  FUNCTION pid_hash_hash(this, hash_id) RESULT(hash)

    CLASS(particle_id_hash), INTENT(IN) :: this
    INTEGER(idkind), INTENT(IN) :: hash_id
    INTEGER(i8) :: hash

    ! Knuth's multiplicative method (sort of)
    hash = INT(MODULO(hash_id * INT(this%hash_gr, idkind), &
        SIZE(this%buckets, KIND=idkind)), idkind) + 1_idkind

  END FUNCTION pid_hash_hash



  !> Test if this hash table holds a given id

  FUNCTION pid_hash_holds(this, test_id) RESULT (holds)

    CLASS(particle_id_hash), INTENT(IN) :: this
    INTEGER(idkind), INTENT(IN) :: test_id
    LOGICAL :: holds
    INTEGER(i8) :: bucket

    IF (.NOT. ALLOCATED(this%buckets)) THEN
      holds = .FALSE.
      RETURN
    END IF

    bucket = this%hash(test_id)
    holds = this%buckets(bucket)%holds(test_id)

  END FUNCTION pid_hash_holds



  !> Add an ID to the hash table

  SUBROUTINE pid_hash_add(this, add_id)

    CLASS(particle_id_hash), INTENT(INOUT) :: this
    INTEGER(idkind), INTENT(IN) :: add_id
    INTEGER(i8) :: bucket

    IF (.NOT. ALLOCATED(this%buckets)) RETURN
    bucket = this%hash(add_id)
    CALL this%buckets(bucket)%add(add_id)
    this%count = this%count + 1

  END SUBROUTINE pid_hash_add



  !> Add those IDs from a list of IDs that correspond to particles on this
  !> processor

  SUBROUTINE pid_hash_add_if_local(this, id_list, sorted)

    CLASS(particle_id_hash), INTENT(INOUT) :: this
    INTEGER(idkind), DIMENSION(:), INTENT(INOUT) :: id_list
    LOGICAL, INTENT(IN), OPTIONAL :: sorted
#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    TYPE(particle), POINTER :: current
    INTEGER :: ispecies

    IF (.NOT. PRESENT(sorted)) THEN
      CALL sort_list(id_list)
    ELSE IF (.NOT. sorted) THEN
      CALL sort_list(id_list)
    END IF
    DO ispecies = 1, n_species
      current => species_list(ispecies)%attached_list%head
      DO WHILE(ASSOCIATED(current))
        IF (in_list(id_list, current%id)) THEN
          CALL this%add(current%id)
          this%count = this%count + 1
        END IF
        current => current%next
      END DO
    END DO
#endif

  END SUBROUTINE pid_hash_add_if_local



  !> Function to add a list of hashes based on a binary file

  SUBROUTINE pid_hash_add_from_file(this, id_file, sorted)

    CLASS(particle_id_hash), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: id_file
    LOGICAL, INTENT(IN), OPTIONAL :: sorted
    INTEGER(idkind), DIMENSION(:), ALLOCATABLE :: id_list
#if defined(PARTICLE_ID)
    INTEGER(KIND = MPI_OFFSET_KIND), PARAMETER :: id_length = 8
    INTEGER, PARAMETER :: mpi_type = MPI_INTEGER8
#else
    INTEGER(KIND = MPI_OFFSET_KIND), PARAMETER :: id_length = 4
    INTEGER, PARAMETER :: mpi_type = MPI_INTEGER4
#endif
    INTEGER(KIND = MPI_OFFSET_KIND) :: nels, nels_remaining
    INTEGER :: handle, errcode, to_read

    ALLOCATE(id_list(this%id_chunk_size))
    IF (rank == 0) PRINT*, 'Opening persistent ID file ', TRIM(id_file)

    CALL MPI_FILE_OPEN(MPI_COMM_SELF, id_file, MPI_MODE_RDONLY, &
        MPI_INFO_NULL, handle, errcode)

    IF (errcode /= 0) THEN
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Persistent ID file ', TRIM(id_file), ' does not exist.'
      END IF
      CALL abort_code(c_err_bad_value)
      RETURN
    END IF

    CALL MPI_FILE_GET_SIZE(handle, nels, errcode)
    nels = nels / id_length
    nels_remaining = nels
    DO WHILE(nels_remaining > 0)
      to_read = INT(MIN(nels_remaining, &
          INT(this%id_chunk_size, MPI_OFFSET_KIND)), i4)
      CALL MPI_FILE_READ(handle, id_list, to_read, mpi_type, &
          MPI_STATUS_IGNORE, errcode)
      CALL this%add_if_local(id_list(1:to_read), sorted)
      nels_remaining = nels_remaining - to_read
    END DO

    CALL MPI_FILE_CLOSE(handle, errcode)
    IF (rank == 0) PRINT*, 'Persistent ID file ', TRIM(id_file), ' read OK'

  END SUBROUTINE pid_hash_add_from_file



  !> Remove an ID from the hash table

  FUNCTION pid_hash_delete(this, del_id) RESULT(holds)

    CLASS(particle_id_hash), INTENT(INOUT) :: this
    INTEGER(idkind), INTENT(IN) :: del_id
    INTEGER(i8) :: bucket
    LOGICAL :: holds

    IF (.NOT. ALLOCATED(this%buckets)) THEN
      holds = .FALSE.
      RETURN
    END IF

    bucket = this%hash(del_id)
    holds = this%buckets(bucket)%delete(del_id)
    IF(holds) this%count = this%count - 1

  END FUNCTION pid_hash_delete



  !> Subroutine to initialise the hash table

  SUBROUTINE pid_hash_init_i8(this, bucket_count, realloc)

    CLASS(particle_id_hash), INTENT(INOUT) :: this
    INTEGER(i8), INTENT(IN) :: bucket_count
    LOGICAL, INTENT(IN), OPTIONAL :: realloc
    TYPE(particle_id_inner_list), DIMENSION(:), ALLOCATABLE :: buckets_old
    INTEGER :: local_count, ierr, ibuck, ipart, seed
    LOGICAL :: should_realloc

    should_realloc = .FALSE.
    IF (PRESENT(realloc)) should_realloc = realloc

    ! Since this RNG is only used for partitioning for the quicksort
    ! algorithm (which it ultimately deterministic) seed from clock always
    IF (.NOT. random_state_set) THEN
      CALL SYSTEM_CLOCK(seed)
      seed = seed + rank
      CALL random_init(seed, random_state)
    END IF

    IF (.NOT. ALLOCATED(this%buckets)) THEN
      CALL MPI_ALLREDUCE(bucket_count, local_count, 1, MPI_INTEGER, MPI_MAX, &
          comm, ierr)

      this%hash_gr = &
          INT(0.5_num * (SQRT(5.0_num) - 1.0_num) * REAL(local_count, num), i8)
      local_count = MAX(local_count, 1000)
      this%count = 0
      this%id_chunk_size = id_global_chunk_size

      ALLOCATE(this%buckets(local_count))
    ELSE
      IF (.NOT. should_realloc) RETURN

      CALL MPI_ALLREDUCE(bucket_count, local_count, 1, MPI_INTEGER, MPI_MAX, &
          comm, ierr)

      local_count = MAX(local_count, 1000)
      IF (local_count == SIZE(this%buckets)) RETURN

      ALLOCATE(buckets_old(SIZE(this%buckets)), SOURCE = this%buckets)
      DEALLOCATE(this%buckets)
      ALLOCATE(this%buckets(local_count))

      DO ibuck = 1, SIZE(buckets_old)
        IF (.NOT. ALLOCATED(buckets_old(ibuck)%list)) CYCLE
        DO ipart = 1, SIZE(buckets_old(ibuck)%list)
          CALL this%add(buckets_old(ibuck)%list(ipart))
        END DO
      END DO
      DEALLOCATE(buckets_old)
    END IF

  END SUBROUTINE pid_hash_init_i8



  !> Subroutine to initialise the hash table

  SUBROUTINE pid_hash_init_i4(this, bucket_count, realloc)

    CLASS(particle_id_hash), INTENT(INOUT) :: this
    INTEGER(i4), INTENT(IN) :: bucket_count
    LOGICAL, INTENT(IN), OPTIONAL :: realloc

    CALL this%init(INT(bucket_count,i8), realloc)

  END SUBROUTINE pid_hash_init_i4



  !> Subroutine to optimise the hash table

  SUBROUTINE pid_optimise(this)

    CLASS(particle_id_hash), INTENT(INOUT) :: this

    CALL this%init(this%count, .TRUE.)

  END SUBROUTINE pid_optimise



  !> Delete all inner lists on destruction

  PURE ELEMENTAL SUBROUTINE pid_hash_destructor(this)

    TYPE(particle_id_hash), INTENT(INOUT) :: this

    IF (.NOT. ALLOCATED(this%buckets)) RETURN
    DEALLOCATE(this%buckets)

  END SUBROUTINE pid_hash_destructor



  !> Get the number of stored hashes

  FUNCTION pidr_get_hash_count(this) RESULT(count)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    INTEGER :: count

    count = 0
    IF (.NOT. ALLOCATED(this%list)) RETURN
    count = SIZE(this%list)

  END FUNCTION pidr_get_hash_count



  !> Get a reference by name to a hash object

  FUNCTION pidr_get_hash_by_name(this, name, must_exist) RESULT(hash_ptr)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(IN), OPTIONAL :: must_exist
    TYPE(particle_id_hash), POINTER :: hash_ptr
    LOGICAL :: exist_local

    exist_local = .FALSE.
    IF (PRESENT(must_exist)) exist_local = must_exist

    hash_ptr => this%get_existing_hash(name)
    IF (.NOT. ASSOCIATED(hash_ptr) .AND. .NOT. exist_local) &
        hash_ptr => this%add_new_hash(name)

  END FUNCTION pidr_get_hash_by_name



  !> Add a new hash to the list

  FUNCTION pidr_add_new_hash(this, name) RESULT(hash_ptr)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name
    TYPE(particle_id_hash_holder), DIMENSION(:), ALLOCATABLE :: temp
    INTEGER :: sz
    CLASS(particle_id_hash), POINTER :: hash_ptr
    LOGICAL :: copyback

    sz = 0
    copyback = ALLOCATED(this%list)
    IF (copyback) THEN
      sz = SIZE(this%list)
      IF (sz > 64) THEN
        hash_ptr => NULL()
        RETURN
      END IF
      ALLOCATE(temp(1:sz), SOURCE = this%list)
      DEALLOCATE(this%list)
    END IF

    ALLOCATE(this%list(sz+1))
    IF (copyback) THEN
      this%list(1:sz) = temp
      DEALLOCATE(temp)
    END IF

    ALLOCATE(this%list(sz+1)%contents)
    this%list(sz+1)%contents%name = name(1:MIN(LEN(name), c_max_string_length))
    hash_ptr => this%list(sz+1)%contents

  END FUNCTION pidr_add_new_hash



  !> Get a pointer to an existing hash by name

  FUNCTION pidr_get_existing_hash_by_name(this, name) RESULT(hash_ptr)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name
    CLASS(particle_id_hash), POINTER :: hash_ptr
    INTEGER :: sz, ilist

    hash_ptr => NULL()
    IF (.NOT. ALLOCATED(this%list)) RETURN

    sz = SIZE(this%list)
    DO ilist = 1, sz
      IF (TRIM(this%list(ilist)%contents%name) == TRIM(name)) THEN
        hash_ptr => this%list(ilist)%contents
        RETURN
      END IF
    END DO

  END FUNCTION pidr_get_existing_hash_by_name



  !> Get a pointer to an existing hash by index

  FUNCTION pidr_get_existing_hash_by_index(this, index) RESULT(hash_ptr)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: index
    CLASS(particle_id_hash), POINTER :: hash_ptr
    INTEGER :: sz

    hash_ptr => NULL()
    IF (.NOT. ALLOCATED(this%list)) RETURN

    sz = SIZE(this%list)
    IF (index < 1 .OR. index > sz) RETURN

    hash_ptr => this%list(index)%contents

  END FUNCTION pidr_get_existing_hash_by_index



  !> Delete an ID from all hashes

  SUBROUTINE pidr_delete_all(this, del_id)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    INTEGER(idkind), INTENT(IN) :: del_id
    INTEGER :: ihash, sz
    LOGICAL :: dummy

    IF (.NOT. ALLOCATED(this%list)) RETURN

    sz = SIZE(this%list)
    DO ihash = 1, sz
      dummy = this%list(ihash)%contents%delete(del_id)
    END DO

  END SUBROUTINE pidr_delete_all



  !> Go through all stored hashes and remove the ID from all of them
  !> then return a bitmask showing which hashes the ID was in

  FUNCTION pidr_delete_and_map(this, test_id) RESULT(hashmap)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    INTEGER(idkind), INTENT(IN) :: test_id
    INTEGER(i8) :: hashmap
    INTEGER :: ihash, sz

    hashmap = 0
    IF (.NOT. ALLOCATED(this%list)) RETURN

    sz = SIZE(this%list)
    DO ihash = 1, sz
      hashmap = ISHFT(hashmap, 1_i8)
      IF (this%list(ihash)%contents%delete(test_id)) hashmap = hashmap + 1_i8
    END DO

  END FUNCTION pidr_delete_and_map



  !> Go through all stored hashes and test if the id is in it
  !> then return a bitmask showing which hashes the ID was in

  FUNCTION pidr_map(this, test_id) RESULT(hashmap)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    INTEGER(idkind), INTENT(IN) :: test_id
    INTEGER(i8) :: hashmap
    INTEGER :: ihash, sz

    hashmap = 0
    IF (.NOT. ALLOCATED(this%list)) RETURN

    sz = SIZE(this%list)
    DO ihash = 1, sz
      hashmap = ISHFT(hashmap, 1_i8)
      IF (this%list(ihash)%contents%holds(test_id)) hashmap = hashmap + 1_i8
    END DO

  END FUNCTION pidr_map



  !> Get the hashmap (bitmask of which hashes the specified ID is contained in)

  SUBROUTINE pidr_add_with_map(this, new_id, hashmap)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    INTEGER(idkind), INTENT(IN) :: new_id
    INTEGER(i8), INTENT(IN) :: hashmap
    INTEGER(i8) :: shifthash
    INTEGER :: ihash, sz

    IF (.NOT. ALLOCATED(this%list)) RETURN
    IF (hashmap == 0) RETURN

    sz = SIZE(this%list)
    shifthash = hashmap
    DO ihash = 1, sz
      IF (IAND(shifthash, 1_i8) /= 0_i8) &
          CALL this%list(ihash)%contents%add(new_id)
      shifthash = ISHFT(shifthash, -1_i8)
    END DO

  END SUBROUTINE pidr_add_with_map



  !> Delete all hash tables on destruction

  PURE ELEMENTAL SUBROUTINE pidr_destructor(this)

    TYPE(particle_id_list_registry), INTENT(INOUT) :: this

    IF (.NOT. ALLOCATED(this%list)) RETURN
    DEALLOCATE(this%list)

  END SUBROUTINE pidr_destructor

END MODULE particle_id_hash_mod
