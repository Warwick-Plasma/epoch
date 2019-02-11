MODULE particle_id_hash_mod

  USE constants
  USE shared_data
  USE random_generator
  USE utilities

  IMPLICIT NONE

  INTEGER, PARAMETER :: id_global_chunk_size = 1000000 ! 8MB chunks at i8
  LOGICAL :: random_state_set = .FALSE.
  TYPE(random_state_type), SAVE :: random_state

#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
  INTEGER, PARAMETER :: hkind = idkind
#else
  INTEGER, PARAMETER :: hkind = MPI_ADDRESS_KIND
#endif

  TYPE :: particle_id_inner_list
    INTEGER(hkind), DIMENSION(:), ALLOCATABLE :: list
    CONTAINS
    PROCEDURE :: holds => pid_inner_list_holds
    PROCEDURE :: add => pid_inner_list_add
    PROCEDURE :: delete => pid_inner_list_delete
#ifdef USE_F03
    FINAL :: pid_inner_list_destructor
#endif
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
    PROCEDURE, PUBLIC :: delete => pid_hash_delete
    PROCEDURE :: init_i8 => pid_hash_init_i8
    PROCEDURE :: init_i4 => pid_hash_init_i4
    GENERIC, PUBLIC :: init => init_i8, init_i4
    PROCEDURE, PUBLIC :: optimise => pid_optimise
#ifdef USE_F03
    FINAL :: pid_hash_destructor
#endif
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
    PROCEDURE, PUBLIC :: reset => pidr_reset
#ifdef USE_F03
    FINAL :: pidr_destructor
#endif
  END TYPE particle_id_list_registry

  TYPE(particle_id_list_registry), SAVE :: id_registry

  PRIVATE
  PUBLIC :: id_registry, particle_id_hash

CONTAINS

  !> Test if particle id is in this bin. Linear search

  FUNCTION pid_inner_list_holds(this, test_id, index_out) RESULT(holds)

    CLASS(particle_id_inner_list), INTENT(IN) :: this
    INTEGER(hkind), INTENT(IN) :: test_id
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
    INTEGER(hkind), INTENT(IN) :: add_id
    INTEGER(hkind), DIMENSION(:), ALLOCATABLE :: temp
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
    INTEGER(hkind), INTENT(IN) :: del_id
    INTEGER(hkind), DIMENSION(:), ALLOCATABLE :: temp
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
    INTEGER(hkind), INTENT(IN) :: hash_id
    INTEGER(i8) :: hash

    ! Knuth's multiplicative method (sort of)
    hash = INT(MODULO(hash_id * INT(this%hash_gr, hkind), &
        SIZE(this%buckets, KIND=hkind)), hkind) + 1_hkind

  END FUNCTION pid_hash_hash



  !> Test if this hash table holds a given id

  FUNCTION pid_hash_holds(this, part) RESULT (holds)

    CLASS(particle_id_hash), INTENT(IN) :: this
    TYPE(particle), POINTER, INTENT(IN) :: part
    LOGICAL :: holds
    INTEGER(hkind) :: test_id
    INTEGER(i8) :: bucket

    IF (.NOT. ALLOCATED(this%buckets)) THEN
      holds = .FALSE.
      RETURN
    END IF

#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    test_id = part%id
#else
    CALL MPI_GET_ADDRESS(part, test_id, errcode)
#endif
    bucket = this%hash(test_id)
    holds = this%buckets(bucket)%holds(test_id)

  END FUNCTION pid_hash_holds



  !> Test if this hash table holds a given id

  FUNCTION pid_hash_holds_hkind(this, test_id) RESULT (holds)

    CLASS(particle_id_hash), INTENT(IN) :: this
    INTEGER(hkind), INTENT(IN) :: test_id
    LOGICAL :: holds
    INTEGER(i8) :: bucket

    IF (.NOT. ALLOCATED(this%buckets)) THEN
      holds = .FALSE.
      RETURN
    END IF

    bucket = this%hash(test_id)
    holds = this%buckets(bucket)%holds(test_id)

  END FUNCTION pid_hash_holds_hkind



  !> Add particle to the hash table

  SUBROUTINE pid_hash_add(this, part)

    CLASS(particle_id_hash), INTENT(INOUT) :: this
    TYPE(particle), POINTER, INTENT(IN) :: part
    INTEGER(KIND=hkind) :: new_id

    IF (.NOT. ALLOCATED(this%buckets)) RETURN

#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    new_id = part%id
#else
    CALL MPI_GET_ADDRESS(part, new_id, errcode)
#endif
    CALL pid_hash_add_hkind(this, new_id)

  END SUBROUTINE pid_hash_add



  !> Add an ID to the hash table

  SUBROUTINE pid_hash_add_hkind(this, add_id)

    CLASS(particle_id_hash), INTENT(INOUT) :: this
    INTEGER(hkind), INTENT(IN) :: add_id
    INTEGER(i8) :: bucket

    IF (.NOT. ALLOCATED(this%buckets)) RETURN

    bucket = this%hash(add_id)
    CALL this%buckets(bucket)%add(add_id)
    this%count = this%count + 1

  END SUBROUTINE pid_hash_add_hkind



  !> Remove an ID from the hash table

  FUNCTION pid_hash_delete(this, del_id) RESULT(holds)

    CLASS(particle_id_hash), INTENT(INOUT) :: this
    INTEGER(hkind), INTENT(IN) :: del_id
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
#ifndef USE_F03
      CALL pid_inner_list_destructor(this%buckets)
#endif
      DEALLOCATE(this%buckets)
      ALLOCATE(this%buckets(local_count))

      DO ibuck = 1, SIZE(buckets_old)
        IF (.NOT. ALLOCATED(buckets_old(ibuck)%list)) CYCLE
        DO ipart = 1, SIZE(buckets_old(ibuck)%list)
          CALL pid_hash_add_hkind(this, buckets_old(ibuck)%list(ipart))
        END DO
      END DO
#ifndef USE_F03
      CALL pid_inner_list_destructor(buckets_old)
#endif
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

#ifdef USE_F03
  PURE ELEMENTAL SUBROUTINE pid_hash_destructor(this)

    TYPE(particle_id_hash), INTENT(INOUT) :: this

    IF (.NOT. ALLOCATED(this%buckets)) RETURN
    DEALLOCATE(this%buckets)

  END SUBROUTINE pid_hash_destructor
#endif



  !> Get the number of stored hashes

  FUNCTION pidr_get_hash_count(this, step) RESULT(count)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    INTEGER, OPTIONAL :: step
    INTEGER :: count
    INTEGER, SAVE :: last_step = -1
    INTEGER, SAVE :: last_count = 0

    count = 0
    IF (.NOT. any_persistent_subset) RETURN

    IF (PRESENT(step)) THEN
      IF (step == last_step) THEN
        count = last_count
        RETURN
      END IF
      last_step = step
    END IF

    IF (ALLOCATED(this%list)) THEN
      count = SIZE(this%list)
    ELSE
      count = 0
    END IF

    CALL MPI_ALLREDUCE(MPI_IN_PLACE, count, 1, MPI_INTEGER, MPI_SUM, comm, &
                       errcode)
    last_count = count

  END FUNCTION pidr_get_hash_count



  !> Get a reference by name to a hash object

  FUNCTION pidr_get_hash_by_name(this, name, must_exist) RESULT(hash_ptr)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(IN), OPTIONAL :: must_exist
    CLASS(particle_id_hash), POINTER :: hash_ptr
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

  SUBROUTINE pidr_delete_all(this, part)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    TYPE(particle), POINTER, INTENT(IN) :: part
    INTEGER(hkind) :: del_id
    INTEGER :: ihash, sz
    LOGICAL :: dummy

    IF (.NOT. any_persistent_subset) RETURN
    IF (.NOT. ALLOCATED(this%list)) RETURN

#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    del_id = part%id
#else
    CALL MPI_GET_ADDRESS(part, del_id, errcode)
#endif
    sz = SIZE(this%list)
    DO ihash = 1, sz
      dummy = this%list(ihash)%contents%delete(del_id)
    END DO

  END SUBROUTINE pidr_delete_all



  !> Go through all stored hashes and remove the ID from all of them
  !> then return a bitmask showing which hashes the ID was in

  FUNCTION pidr_delete_and_map(this, test_id) RESULT(hashmap)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    INTEGER(hkind), INTENT(IN) :: test_id
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

  FUNCTION pidr_map(this, part) RESULT(hashmap)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    TYPE(particle), POINTER, INTENT(IN) :: part
    INTEGER(i8) :: hashmap
    INTEGER(hkind) :: test_id
    INTEGER :: ihash, sz

    hashmap = 0
    IF (.NOT. ALLOCATED(this%list)) RETURN

#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    test_id = part%id
#else
    CALL MPI_GET_ADDRESS(part, test_id, errcode)
#endif

    sz = SIZE(this%list)
    DO ihash = 1, sz
      hashmap = ISHFT(hashmap, 1_i8)
      IF (pid_hash_holds_hkind(this%list(ihash)%contents, test_id)) &
          hashmap = hashmap + 1_i8
    END DO

  END FUNCTION pidr_map



  !> Get the hashmap (bitmask of which hashes the specified ID is contained in)

  SUBROUTINE pidr_add_with_map(this, part, hashmap)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this
    TYPE(particle), POINTER, INTENT(IN) :: part
    INTEGER(i8), INTENT(IN) :: hashmap
    INTEGER(hkind) :: new_id
    INTEGER(i8) :: shifthash
    INTEGER :: ihash, sz

    IF (.NOT. ALLOCATED(this%list)) RETURN
    IF (hashmap == 0) RETURN

#if defined(PARTICLE_ID) || defined(PARTICLE_ID4)
    new_id = part%id
#else
    CALL MPI_GET_ADDRESS(part, new_id, errcode)
#endif

    sz = SIZE(this%list)
    shifthash = hashmap
    DO ihash = 1, sz
      IF (IAND(shifthash, 1_i8) /= 0_i8) &
          CALL pid_hash_add_hkind(this%list(ihash)%contents, new_id)
      shifthash = ISHFT(shifthash, -1_i8)
    END DO

  END SUBROUTINE pidr_add_with_map



  !> Delete all hash tables

  PURE ELEMENTAL SUBROUTINE pidr_reset(this)

    CLASS(particle_id_list_registry), INTENT(INOUT) :: this

    IF (.NOT. ALLOCATED(this%list)) RETURN
    DEALLOCATE(this%list)

  END SUBROUTINE pidr_reset



  !> Delete all hash tables on destruction

#ifdef USE_F03
  PURE ELEMENTAL SUBROUTINE pidr_destructor(this)

    TYPE(particle_id_list_registry), INTENT(INOUT) :: this

    CALL pidr_reset(this)

  END SUBROUTINE pidr_destructor
#endif

END MODULE particle_id_hash_mod
