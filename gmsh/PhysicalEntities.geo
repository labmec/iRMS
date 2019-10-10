
Macro DrawBoundaries

If (dimension == 3)



Physical Volume("d_rock") = {v1};
Physical Volume("d_wellbore_p") = {well_p_v_regions[]};
Physical Volume("d_wellbore_i") = {well_i_v_regions[]};
Physical Surface("bc_bottom") = {res_B[]};
Physical Surface("bc_top") = {res_T[]};
Physical Surface("bc_south") = {res_S[]};
Physical Surface("bc_east") = {res_E[]};
Physical Surface("bc_north") = {res_N[]};
Physical Surface("bc_west") = {res_W[]};

// Tagging boundary conditions for prodution wells
Physical Surface("bc_well_p_lids") = well_p_lids[];
Physical Surface("well_i_lids") = well_i_lids[];
Physical Surface("bc_producers") = well_p_bores[];
Physical Surface("bc_injectors") = well_i_bores[];


// side-burden
If (geomechanicQ == 1)
Physical Volume("side_burden") = {v2};
Physical Surface("side_burden_bottom") = {sb_B[]};
Physical Surface("side_burden_top") = {sb_T[]};
Physical Surface("side_burden_South") = {sb_S[]};
Physical Surface("side_burden_East") = {sb_E[]};
Physical Surface("side_burden_North") = {sb_N[]};
Physical Surface("side_burden_West") = {sb_W[]};
EndIf



Else



Physical Surface("d_rock") = {s1};
Physical Surface("d_Wellbore_p") = {well_p_v_regions[]};
Physical Surface("d_Wellbore_i") = {well_i_v_regions[]};
Physical Line("bc_south") = {res_S[]};
Physical Line("bc_West") = {res_W[]};
Physical Line("bc_north") = {res_N[]};
Physical Line("bc_East") = {res_E[]};

// Tagging boundary conditions for prodution wells
Physical Line("bc_well_p_lids") = well_p_lids[];
Physical Line("bc_well_i_lids") = well_i_lids[];
Physical Line("bc_producers") = well_p_bores[];
Physical Line("bc_injectors") = well_i_bores[];

// side-burden
If (geomechanicQ == 1)
Physical Surface("side_burden") = {s2};
Physical Line("side_burden_south") = {sb_S[]};
Physical Line("side_burden_West") = {sb_W[]};
Physical Line("side_burden_north") = {sb_N[]};
Physical Line("side_burden_East") = {sb_E[]};
EndIf

EndIf

Return

