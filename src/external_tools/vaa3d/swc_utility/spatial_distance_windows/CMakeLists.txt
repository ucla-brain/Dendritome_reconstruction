message("external tool: spatial distance (from Vaa3D)")
set(VAA3D_SPATIAL_DISTANCE_SRC
    xyz.hpp
    neuron_tree.hpp
    swc_file.hpp
    v_neuronswc.hpp
    neuron_sim_scores.hpp
    swc_file.cpp
    v_neuronswc.cpp
    neuron_sim_scores.cpp
    spatial_distance.cpp
    )

add_executable(spatial_distance ${VAA3D_SPATIAL_DISTANCE_SRC})

