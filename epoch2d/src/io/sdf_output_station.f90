MODULE sdf_output_station

  USE sdf_output_station_r4
  USE sdf_output_station_r8

  IMPLICIT NONE

  INTERFACE sdf_write_station_header
    MODULE PROCEDURE &
        write_station_header_1d_r4, &
        write_station_header_2d_r4, &
        write_station_header_3d_r4, &
        write_station_header_1d_r8, &
        write_station_header_2d_r8, &
        write_station_header_3d_r8
  END INTERFACE sdf_write_station_header

  INTERFACE sdf_write_station_array
    MODULE PROCEDURE &
        write_station_array_r4_r4, &
        write_station_array_r8_r4, &
        write_station_array_r4_r8, &
        write_station_array_r8_r8
  END INTERFACE sdf_write_station_array

END MODULE sdf_output_station
