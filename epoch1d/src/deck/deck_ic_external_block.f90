MODULE deck_ic_external_block

  USE shared_data
  USE strings_advanced
  USE simple_io

  IMPLICIT NONE

  SAVE

  INTEGER(KIND=MPI_OFFSET_KIND) :: offset
CONTAINS

  SUBROUTINE start_external
    offset=0
  END SUBROUTINE start_external

  FUNCTION handle_ic_external_species_deck(species_id,element,value)
    CHARACTER(*),INTENT(IN) :: element,value
    INTEGER,INTENT(IN) :: species_id
    INTEGER :: handle_ic_external_species_deck
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset =0

    handle_ic_external_species_deck=c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element,"offset")) THEN
       offset=as_long_integer_simple(value,handle_ic_external_species_deck)
       RETURN
    ENDIF

    IF (species_id .LT. 0 .OR. species_id .GT. n_species) THEN
       IF (rank .EQ. 0) PRINT *,"Attempting to set non-existent species initial conditions. Ignoring."
       RETURN
    ENDIF

    IF (str_cmp(element,"minrho")) THEN
       initial_conditions(species_id)%minrho=as_real(value,handle_ic_external_species_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"maxrho")) THEN
       initial_conditions(species_id)%maxrho=as_real(value,handle_ic_external_species_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"rho")) THEN
       CALL load_single_array_from_data_file(value,initial_conditions(species_id)%rho,offset,handle_ic_external_species_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"temp")) THEN
       CALL load_single_array_from_data_file(value,initial_conditions(species_id)%temp(:,1),offset&
            ,handle_ic_external_species_deck)
       initial_conditions(species_id)%temp(-1:nx+2,2)=initial_conditions(species_id)%temp(-1:nx+2,1)
       initial_conditions(species_id)%temp(-1:nx+2,3)=initial_conditions(species_id)%temp(-1:nx+2,1)
       RETURN
    ENDIF


    IF (str_cmp(element,"temp_x")) THEN
       CALL load_single_array_from_data_file(value,initial_conditions(species_id)%temp(:,1),offset&
            ,handle_ic_external_species_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"temp_y")) THEN
       CALL load_single_array_from_data_file(value,initial_conditions(species_id)%temp(:,2),offset&
            ,handle_ic_external_species_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"temp_z")) THEN
       CALL load_single_array_from_data_file(value,initial_conditions(species_id)%temp(:,3),offset&
            ,handle_ic_external_species_deck)
       RETURN
    ENDIF
    handle_ic_external_species_deck=c_err_unknown_element

  END FUNCTION handle_ic_external_species_deck

  FUNCTION handle_ic_external_fields_deck(element,value)
    CHARACTER(*),INTENT(IN) :: element,value
    INTEGER :: handle_ic_external_fields_deck
    INTEGER(KIND=MPI_OFFSET_KIND) :: offset =0

    handle_ic_external_fields_deck=c_err_none
    IF (element .EQ. blank .OR. value .EQ. blank) RETURN

    IF (str_cmp(element,"ex")) THEN
       CALL load_single_array_from_data_file(value,ex,offset,handle_ic_external_fields_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"ey")) THEN
       CALL load_single_array_from_data_file(value,ey,offset,handle_ic_external_fields_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"ez")) THEN
       CALL load_single_array_from_data_file(value,ez,offset,handle_ic_external_fields_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"bx")) THEN
       CALL load_single_array_from_data_file(value,bx,offset,handle_ic_external_fields_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"by")) THEN
       CALL load_single_array_from_data_file(value,by,offset,handle_ic_external_fields_deck)
       RETURN
    ENDIF

    IF (str_cmp(element,"bz")) THEN
       CALL load_single_array_from_data_file(value,bz,offset,handle_ic_external_fields_deck)
       RETURN
    ENDIF

  END FUNCTION handle_ic_external_fields_deck

END MODULE deck_ic_external_block
