<map version="0.9.0">
<!-- To view this file, download free mind mapping software FreeMind from http://freemind.sourceforge.net -->
<node CREATED="1307791262350" ID="ID_319356386" MODIFIED="1319668026996" TEXT="main (main.cpp)">
<node CREATED="1307793288498" ID="ID_1385423889" MODIFIED="1319667922897" POSITION="right" TEXT="read_defines (main.cpp)"/>
<node CREATED="1307791279195" ID="ID_751379187" MODIFIED="1319667990672" POSITION="right" TEXT="initialization  (main.cpp)">
<node CREATED="1307791551494" ID="ID_1793305870" MODIFIED="1320495217627" TEXT="communication_initialization (mpi.cpp)"/>
<node CREATED="1307791564589" ID="ID_1926601743" MODIFIED="1319668182517" TEXT="global_to_local (main.cpp)"/>
<node CREATED="1307791573275" ID="ID_664777153" MODIFIED="1320495525507" TEXT="device_initialization (cpu/gpu)">
<node CREATED="1307803548226" ID="ID_733027286" MODIFIED="1319668247329" TEXT="&#x438;&#x43d;&#x438;&#x446;&#x438;&#x430;&#x43b;&#x438;&#x437;&#x430;&#x446;&#x438;&#x44f; GPU (gpu)"/>
<node CREATED="1307803566459" ID="ID_1203476732" MODIFIED="1319668239957" TEXT="&#x43a;&#x43e;&#x43f;&#x438;&#x440;&#x43e;&#x432;&#x430;&#x43d;&#x438;&#x435; &#x43a;&#x43e;&#x43d;&#x441;&#x442;&#x430;&#x43d;&#x442; &#x432; &#x43f;&#x430;&#x43c;&#x44f;&#x442;&#x44c; GPU (gpu)"/>
</node>
<node CREATED="1307791582136" ID="ID_1086307714" MODIFIED="1319668256113" TEXT="memory_allocation  (main.cpp)">
<node CREATED="1307792209741" ID="ID_81072915" MODIFIED="1319668194151" TEXT="host_memory_allocation (main.cpp)"/>
<node CREATED="1307792218142" ID="ID_1763565157" MODIFIED="1319668210394" TEXT="device_memory_allocation (cpu/gpu)"/>
</node>
<node CREATED="1307791592497" ID="ID_315198654" MODIFIED="1319668309098" TEXT="restore (main.cpp)">
<node CREATED="1307792296114" ID="ID_606683520" MODIFIED="1319668312743" TEXT="&#x447;&#x442;&#x435;&#x43d;&#x438;&#x435; &#x438;&#x437; &#x444;&#x430;&#x439;&#x43b;&#x430; (P1, S2, x, y, z, roS1_old, roS2_old, media) (main.cpp)"/>
</node>
<node CREATED="1307791615149" ID="ID_792452932" MODIFIED="1319668094652" TEXT="data_initialization (main.cpp)">
<node CREATED="1307792367068" ID="ID_1436541564" MODIFIED="1319668109056" TEXT="&#x43f;&#x440;&#x438;&#x441;&#x432;&#x43e;&#x435;&#x43d;&#x438;&#x435; (S2, P1, x, y, z, media)  (main.cpp)"/>
<node CREATED="1307792392220" ID="ID_316773299" MODIFIED="1319668111261" TEXT="&#x437;&#x430;&#x434;&#x430;&#x43d;&#x438;&#x435; &#x441;&#x442;&#x440;&#x443;&#x43a;&#x442;&#x443;&#x440;&#x44b; &#x440;&#x435;&#x437;&#x435;&#x440;&#x432;&#x443;&#x430;&#x440;&#x430; (main.cpp)"/>
</node>
<node CREATED="1307791626633" ID="ID_192250718" MODIFIED="1319668337988" TEXT="load_data_to_device (P1, S2, roS1_old, roS2_old, media) (cpu/gpu)"/>
</node>
<node CREATED="1307791370055" ID="ID_1982838687" MODIFIED="1319668002014" POSITION="right" TEXT="time_step_function (main.cpp)">
<node CREATED="1307792478624" ID="ID_1977192121" MODIFIED="1320495499406" TEXT="P_S_exchange (mpi.cpp)">
<node CREATED="1307792682077" ID="ID_1596445838" MODIFIED="1320495676330" TEXT="exchange (P_w, S_n) (mpi.cpp)">
<node CREATED="1307792741446" ID="ID_1468306152" MODIFIED="1320495581904" TEXT="load_data_to_host (P_w) (cpu/gpu)"/>
<node CREATED="1307792791859" ID="ID_1097133708" MODIFIED="1320495590745" TEXT="right_send_recv (mpi.cpp)"/>
<node CREATED="1307792807911" ID="ID_515651811" MODIFIED="1320495599522" TEXT="left_recv_send (mpi.cpp)"/>
<node CREATED="1307792829731" ID="ID_1266197297" MODIFIED="1320495609573" TEXT="load_data_to_device(P_w) (cpu/gpu)"/>
</node>
</node>
<node CREATED="1307792499708" ID="ID_145217363" MODIFIED="1321218204248" TEXT="ro_P_Xi_calculation (cpu/gpu)">
<node CREATED="1321218186398" ID="ID_190010606" MODIFIED="1321218219648" TEXT="ro_P_Xi_calculation_kernel(cpu-project/gpu)"/>
</node>
<node CREATED="1307792513407" ID="ID_216870866" MODIFIED="1320495505135" TEXT="P_ro_Xi_exchange (mpi.cpp)">
<node CREATED="1307792633991" ID="ID_295742772" MODIFIED="1320495412084" TEXT="exchange (P_n, ro_w, ro_n, Xi_w, Xi_n)"/>
</node>
<node CREATED="1307792525392" ID="ID_1124813596" MODIFIED="1320495538930" TEXT="u_calculation (cpu/gpu)">
<node CREATED="1307793160114" ID="ID_1407139979" MODIFIED="1307793168601" TEXT="assign_u_kernel"/>
</node>
<node CREATED="1307792533432" ID="ID_1979749368" MODIFIED="1320495507654" TEXT="u_exchange (mpi.cpp)">
<node CREATED="1307792657995" ID="ID_422162669" MODIFIED="1320495470504" TEXT="exchange (ux_w, ux_n)"/>
</node>
<node CREATED="1307792538878" ID="ID_1310497056" MODIFIED="1320495540977" TEXT="roS_calculation (cpu/gpu)">
<node CREATED="1307793125965" ID="ID_1070652216" MODIFIED="1307793138518" TEXT="assign_rS_kernel"/>
<node CREATED="1307793139269" ID="ID_998884838" MODIFIED="1307793153218" TEXT="assign_rS_kernel_nr"/>
</node>
<node CREATED="1307792547022" ID="ID_686920280" MODIFIED="1321218226487" TEXT="P_S_calculation (cpu/gpu)">
<node CREATED="1307793174456" ID="ID_1077397314" MODIFIED="1321218245234" TEXT="Newton_method_kernel (cpu-project/gpu)"/>
</node>
<node CREATED="1307792555709" ID="ID_1157062549" MODIFIED="1321218249388" TEXT="boundary_conditions (cpu/gpu)">
<node CREATED="1307793074539" ID="ID_13790351" MODIFIED="1321218267495" TEXT="&#x43f;&#x43e; S_n (cpu-project/gpu)"/>
<node CREATED="1307793082335" ID="ID_247231528" MODIFIED="1321218273348" TEXT="&#x43f;&#x43e; P_w (cpu-project/gpu)"/>
</node>
</node>
<node CREATED="1307791391402" ID="ID_695841967" MODIFIED="1319668006704" POSITION="right" TEXT="save_data_plots (main.cpp)">
<node CREATED="1307791757510" ID="ID_1950853460" MODIFIED="1307791826109" TEXT="load_data_to_host (P1, S2, u2xyz)"/>
<node CREATED="1307791852356" ID="ID_513344022" MODIFIED="1319668454046" TEXT="print_plots_top (main.cpp)"/>
<node CREATED="1307791861461" ID="ID_288962519" MODIFIED="1319668451303" TEXT="print_plots (main.cpp)"/>
</node>
<node CREATED="1307791981139" ID="ID_1622351155" MODIFIED="1319668010341" POSITION="right" TEXT="save (main.cpp)">
<node CREATED="1307791984768" ID="ID_827264247" MODIFIED="1320495630658" TEXT="load_data_to_host (roS1_old, roS2_old) (cpu/gpu)"/>
<node CREATED="1307792029222" ID="ID_1145686532" MODIFIED="1320495656723" TEXT="&#x437;&#x430;&#x43f;&#x438;&#x441;&#x44c; &#x432; &#x444;&#x430;&#x439;&#x43b; (P_w, S_n, x, y, z, roS_w_old, roS_n_old, media) (main.cpp)"/>
</node>
<node CREATED="1307791444495" ID="ID_1398112335" MODIFIED="1319668021076" POSITION="right" TEXT="finalization (main.cpp)">
<node CREATED="1307791685505" ID="ID_1093742610" MODIFIED="1319668390375" TEXT="memory_free (main.cpp)">
<node CREATED="1307791724431" ID="ID_1791646214" MODIFIED="1319668413665" TEXT="host_memory_free (main.cpp)"/>
<node CREATED="1307791733588" ID="ID_692188628" MODIFIED="1319668434942" TEXT="device_memory_free (cpu/gpu)"/>
</node>
<node CREATED="1307791697237" ID="ID_671746965" MODIFIED="1319668409000" TEXT="communication_finalization (mpi.cpp)"/>
<node CREATED="1321045894405" ID="ID_271479106" MODIFIED="1321045920898" TEXT="device_finalization (gpu/cpu)"/>
</node>
</node>
</map>
