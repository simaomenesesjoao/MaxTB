/*
 * rhombic.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Oct 6 2022, 18:05 by COMSOL 5.5.0.359. */
public class comsol_freq1p8L20 {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/workspace/sjoao");

    model.label("rhombic.mph");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 3);

    model.result().table().create("evl3", "Table");

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").lengthUnit("nm");
    model.component("comp1").geom("geom1").create("wp1", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp1").label("xy1");
    model.component("comp1").geom("geom1").feature("wp1").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp1").set("normalvector", new int[]{1, 1, 0});
    model.component("comp1").geom("geom1").feature("wp1").set("normalcoord", new int[]{1, 1, 0});
    model.component("comp1").geom("geom1").feature("wp1").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp1").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp1").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp1").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp2", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp2").label("xy2");
    model.component("comp1").geom("geom1").feature("wp2").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp2").set("normalvector", new int[]{1, -1, 0});
    model.component("comp1").geom("geom1").feature("wp2").set("normalcoord", new int[]{1, -1, 0});
    model.component("comp1").geom("geom1").feature("wp2").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp2").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp2").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp2").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp3", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp3").label("xy3");
    model.component("comp1").geom("geom1").feature("wp3").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp3").set("normalvector", new int[]{-1, -1, 0});
    model.component("comp1").geom("geom1").feature("wp3").set("normalcoord", new int[]{-1, -1, 0});
    model.component("comp1").geom("geom1").feature("wp3").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp3").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp3").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp3").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp4", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp4").label("xy4");
    model.component("comp1").geom("geom1").feature("wp4").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp4").set("normalvector", new int[]{-1, 1, 0});
    model.component("comp1").geom("geom1").feature("wp4").set("normalcoord", new int[]{-1, 1, 0});
    model.component("comp1").geom("geom1").feature("wp4").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp4").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp4").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp4").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp5", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp5").label("xz1");
    model.component("comp1").geom("geom1").feature("wp5").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp5").set("normalvector", new int[]{1, 0, 1});
    model.component("comp1").geom("geom1").feature("wp5").set("normalcoord", new int[]{1, 1, 1});
    model.component("comp1").geom("geom1").feature("wp5").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp5").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp5").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp5").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp6", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp6").label("xz2");
    model.component("comp1").geom("geom1").feature("wp6").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp6").set("normalvector", new int[]{1, 0, -1});
    model.component("comp1").geom("geom1").feature("wp6").set("normalcoord", new int[]{1, 0, -1});
    model.component("comp1").geom("geom1").feature("wp6").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp6").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp6").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp6").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp7", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp7").label("xz3");
    model.component("comp1").geom("geom1").feature("wp7").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp7").set("normalvector", new int[]{-1, 0, -1});
    model.component("comp1").geom("geom1").feature("wp7").set("normalcoord", new int[]{-1, 0, -1});
    model.component("comp1").geom("geom1").feature("wp7").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp7").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp7").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp7").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp8", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp8").label("xz4");
    model.component("comp1").geom("geom1").feature("wp8").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp8").set("normalvector", new int[]{-1, 0, 1});
    model.component("comp1").geom("geom1").feature("wp8").set("normalcoord", new int[]{-1, 0, 1});
    model.component("comp1").geom("geom1").feature("wp8").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp8").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp8").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp8").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp9", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp9").label("yz1");
    model.component("comp1").geom("geom1").feature("wp9").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp9").set("normalvector", new int[]{0, 1, 1});
    model.component("comp1").geom("geom1").feature("wp9").set("normalcoord", new int[]{0, 1, 1});
    model.component("comp1").geom("geom1").feature("wp9").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp9").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp9").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp9").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp10", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp10").label("yz2");
    model.component("comp1").geom("geom1").feature("wp10").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp10").set("normalvector", new int[]{0, 1, -1});
    model.component("comp1").geom("geom1").feature("wp10").set("normalcoord", new int[]{0, 1, -1});
    model.component("comp1").geom("geom1").feature("wp10").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp10").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp10").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp10").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp11", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp11").label("yz3");
    model.component("comp1").geom("geom1").feature("wp11").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp11").set("normalvector", new int[]{0, -1, -1});
    model.component("comp1").geom("geom1").feature("wp11").set("normalcoord", new int[]{0, -1, -1});
    model.component("comp1").geom("geom1").feature("wp11").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp11").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp11").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp11").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("wp12", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp12").label("yz4");
    model.component("comp1").geom("geom1").feature("wp12").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp12").set("normalvector", new int[]{0, -1, 1});
    model.component("comp1").geom("geom1").feature("wp12").set("normalcoord", new int[]{0, -1, 1});
    model.component("comp1").geom("geom1").feature("wp12").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp12").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp12").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp12").geom().feature("sq1").set("size", 6);
    model.component("comp1").geom("geom1").create("ext1", "Extrude");
    model.component("comp1").geom("geom1").feature("ext1").set("workplane", "wp1");
    model.component("comp1").geom("geom1").feature("ext1").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext1").selection("input").set("wp1");
    model.component("comp1").geom("geom1").create("ext2", "Extrude");
    model.component("comp1").geom("geom1").feature("ext2").set("workplane", "wp2");
    model.component("comp1").geom("geom1").feature("ext2").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext2").selection("input").set("wp2");
    model.component("comp1").geom("geom1").create("ext3", "Extrude");
    model.component("comp1").geom("geom1").feature("ext3").set("workplane", "wp3");
    model.component("comp1").geom("geom1").feature("ext3").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext3").selection("input").set("wp3");
    model.component("comp1").geom("geom1").create("ext4", "Extrude");
    model.component("comp1").geom("geom1").feature("ext4").set("workplane", "wp4");
    model.component("comp1").geom("geom1").feature("ext4").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext4").selection("input").set("wp4");
    model.component("comp1").geom("geom1").create("ext5", "Extrude");
    model.component("comp1").geom("geom1").feature("ext5").set("workplane", "wp5");
    model.component("comp1").geom("geom1").feature("ext5").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext5").selection("input").set("wp5");
    model.component("comp1").geom("geom1").create("ext6", "Extrude");
    model.component("comp1").geom("geom1").feature("ext6").set("workplane", "wp6");
    model.component("comp1").geom("geom1").feature("ext6").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext6").selection("input").set("wp6");
    model.component("comp1").geom("geom1").create("ext7", "Extrude");
    model.component("comp1").geom("geom1").feature("ext7").set("workplane", "wp7");
    model.component("comp1").geom("geom1").feature("ext7").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext7").selection("input").set("wp7");
    model.component("comp1").geom("geom1").create("ext8", "Extrude");
    model.component("comp1").geom("geom1").feature("ext8").set("workplane", "wp8");
    model.component("comp1").geom("geom1").feature("ext8").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext8").selection("input").set("wp8");
    model.component("comp1").geom("geom1").create("ext9", "Extrude");
    model.component("comp1").geom("geom1").feature("ext9").set("workplane", "wp9");
    model.component("comp1").geom("geom1").feature("ext9").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext9").selection("input").set("wp9");
    model.component("comp1").geom("geom1").create("ext10", "Extrude");
    model.component("comp1").geom("geom1").feature("ext10").set("workplane", "wp10");
    model.component("comp1").geom("geom1").feature("ext10").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext10").selection("input").set("wp10");
    model.component("comp1").geom("geom1").create("ext11", "Extrude");
    model.component("comp1").geom("geom1").feature("ext11").set("workplane", "wp11");
    model.component("comp1").geom("geom1").feature("ext11").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext11").selection("input").set("wp11");
    model.component("comp1").geom("geom1").create("ext12", "Extrude");
    model.component("comp1").geom("geom1").feature("ext12").setIndex("distance", "-5", 0);
    model.component("comp1").geom("geom1").feature("ext12").selection("input").set("wp12");
    model.component("comp1").geom("geom1").create("int1", "Intersection");
    model.component("comp1").geom("geom1").feature("int1").selection("input")
         .set("ext1", "ext10", "ext11", "ext12", "ext2", "ext3", "ext4", "ext5", "ext6", "ext7", "ext8", "ext9");
    model.component("comp1").geom("geom1").create("sca1", "Scale");
    model.component("comp1").geom("geom1").feature("sca1").set("factor", 0.75);
    model.component("comp1").geom("geom1").feature("sca1").selection("input").set("int1");
    model.component("comp1").geom("geom1").create("blk1", "Block");
    model.component("comp1").geom("geom1").feature("blk1").set("pos", new int[]{0, 0, 0});
    model.component("comp1").geom("geom1").feature("blk1").set("base", "center");
    model.component("comp1").geom("geom1").feature("blk1").set("size", new int[]{20, 20, 20});
    model.component("comp1").geom("geom1").run();

    model.component("comp1").selection().create("sel1", "Explicit");
    model.component("comp1").selection().create("sel2", "Explicit");

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material().create("mat2", "Common");
    model.component("comp1").material("mat1").selection().set(2);
    model.component("comp1").material("mat2").selection().set(1);

    model.component("comp1").physics().create("es", "Electrostatics", "geom1");
    model.component("comp1").physics("es").create("pot1", "ElectricPotential", 2);
    model.component("comp1").physics("es").feature("pot1").selection().set(1, 2, 3, 4, 5, 18);

    model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");

    model.result().table("evl3").label("Evaluation 3D");
    model.result().table("evl3").comments("Interactive 3D values");

    model.component("comp1").view("view1").set("transparency", true);
    model.component("comp1").view("view2").axis().set("xmin", -9.087827682495117);
    model.component("comp1").view("view2").axis().set("xmax", 15.439041137695312);
    model.component("comp1").view("view2").axis().set("ymin", -11.766670227050781);
    model.component("comp1").view("view2").axis().set("ymax", 6.3875651359558105);

    model.component("comp1").material("mat1").propertyGroup("def")
         .set("relpermittivity", new String[]{"-14.584299999999999+0.6876i", "0", "0", "0", "-14.584299999999999+0.6876i", "0", "0", "0", "-14.584299999999999+0.6876i"});
    model.component("comp1").material("mat1").propertyGroup("def").set("relpermittivity_symmetry", "0");
    model.component("comp1").material("mat2").propertyGroup("def")
         .set("relpermittivity", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat2").propertyGroup("def").set("relpermittivity_symmetry", "0");

    model.component("comp1").physics("es").feature("pot1").set("V0", "z");

    model.component("comp1").mesh("mesh1").feature("size").set("hauto", 1);
    model.component("comp1").mesh("mesh1").run();

    model.study().create("std1");
    model.study("std1").create("stat", "Stationary");

    model.sol().create("sol1");
    model.sol("sol1").study("std1");
    model.sol("sol1").attach("std1");
    model.sol("sol1").create("st1", "StudyStep");
    model.sol("sol1").create("v1", "Variables");
    model.sol("sol1").create("s1", "Stationary");
    model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
    model.sol("sol1").feature("s1").create("i1", "Iterative");
    model.sol("sol1").feature("s1").feature("i1").create("mg1", "Multigrid");
    model.sol("sol1").feature("s1").feature().remove("fcDef");

    model.result().create("pg1", "PlotGroup3D");
    model.result("pg1").create("mslc1", "Multislice");
    model.result().export().create("data1", "Data");

    model.sol("sol1").attach("std1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("prefun", "amg");
    model.sol("sol1").runAll();

    model.result("pg1").label("Electric Potential (es)");
    model.result("pg1").set("frametype", "spatial");
    model.result("pg1").feature("mslc1").set("colortable", "RainbowLight");
    model.result("pg1").feature("mslc1").set("resolution", "normal");
    model.result().export("data1").set("expr", new String[]{"V"});
    model.result().export("data1").set("unit", new String[]{"V"});
    model.result().export("data1").set("descr", new String[]{"Electric potential"});
    model.result().export("data1").set("filename", "/home/sjoao/data/rhombic/pot_freq1p8L20_rad8.0.dat");
    model.result().export("data1").set("location", "file");
    model.result().export("data1").set("coordfilename", "/home/sjoao/data/rhombic/positions_rad8.0.dat");
    model.result().export("data1").run();

    return model;
  }

  public static void main(String[] args) {
    run();
  }

}
