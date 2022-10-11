/*
 * octa2.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Oct 3 2022, 17:19 by COMSOL 5.5.0.359. */
public class comsol_octa_freq2p7 {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/home/sjoao/Desktop");

    model.label("octa2.mph");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 3);

    model.result().table().create("evl3", "Table");

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").lengthUnit("nm");
    model.component("comp1").geom("geom1").create("wp1", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp1").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp1").set("normalvector", new int[]{1, 1, -1});
    model.component("comp1").geom("geom1").feature("wp1").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp1").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp1").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp1").geom().feature("sq1").set("size", 10);
    model.component("comp1").geom("geom1").create("wp2", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp2").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp2").set("normalvector", new int[]{1, -1, -1});
    model.component("comp1").geom("geom1").feature("wp2").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp2").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp2").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp2").geom().feature("sq1").set("size", 10);
    model.component("comp1").geom("geom1").create("wp3", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp3").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp3").set("normalvector", new int[]{-1, 1, -1});
    model.component("comp1").geom("geom1").feature("wp3").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp3").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp3").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp3").geom().feature("sq1").set("size", 10);
    model.component("comp1").geom("geom1").create("wp4", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp4").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp4").set("normalvector", new int[]{-1, -1, -1});
    model.component("comp1").geom("geom1").feature("wp4").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp4").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp4").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp4").geom().feature("sq1").set("size", 10);
    model.component("comp1").geom("geom1").create("ext1", "Extrude");
    model.component("comp1").geom("geom1").feature("ext1").set("workplane", "wp1");
    model.component("comp1").geom("geom1").feature("ext1").setIndex("distance", "-10", 0);
    model.component("comp1").geom("geom1").feature("ext1").selection("input").set("wp1");
    model.component("comp1").geom("geom1").create("ext2", "Extrude");
    model.component("comp1").geom("geom1").feature("ext2").set("workplane", "wp2");
    model.component("comp1").geom("geom1").feature("ext2").setIndex("distance", "-10", 0);
    model.component("comp1").geom("geom1").feature("ext2").selection("input").set("wp2");
    model.component("comp1").geom("geom1").create("ext3", "Extrude");
    model.component("comp1").geom("geom1").feature("ext3").set("workplane", "wp3");
    model.component("comp1").geom("geom1").feature("ext3").setIndex("distance", "-10", 0);
    model.component("comp1").geom("geom1").feature("ext3").selection("input").set("wp3");
    model.component("comp1").geom("geom1").create("ext4", "Extrude");
    model.component("comp1").geom("geom1").feature("ext4").setIndex("distance", "-10", 0);
    model.component("comp1").geom("geom1").feature("ext4").selection("input").set("wp4");
    model.component("comp1").geom("geom1").create("wp5", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp5").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp5").set("normalvector", new int[]{1, 1, -1});
    model.component("comp1").geom("geom1").feature("wp5").set("normalcoord", new int[]{0, 0, 3});
    model.component("comp1").geom("geom1").feature("wp5").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp5").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp5").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp5").geom().feature("sq1").set("size", 10);
    model.component("comp1").geom("geom1").create("wp6", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp6").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp6").set("normalvector", new int[]{1, -1, -1});
    model.component("comp1").geom("geom1").feature("wp6").set("normalcoord", new int[]{0, 0, 3});
    model.component("comp1").geom("geom1").feature("wp6").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp6").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp6").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp6").geom().feature("sq1").set("size", 10);
    model.component("comp1").geom("geom1").create("wp7", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp7").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp7").set("normalvector", new int[]{-1, 1, -1});
    model.component("comp1").geom("geom1").feature("wp7").set("normalcoord", new int[]{0, 0, 3});
    model.component("comp1").geom("geom1").feature("wp7").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp7").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp7").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp7").geom().feature("sq1").set("size", 10);
    model.component("comp1").geom("geom1").create("wp8", "WorkPlane");
    model.component("comp1").geom("geom1").feature("wp8").set("planetype", "normalvector");
    model.component("comp1").geom("geom1").feature("wp8").set("normalvector", new int[]{-1, -1, -1});
    model.component("comp1").geom("geom1").feature("wp8").set("normalcoord", new int[]{0, 0, 3});
    model.component("comp1").geom("geom1").feature("wp8").set("unite", true);
    model.component("comp1").geom("geom1").feature("wp8").geom().create("sq1", "Square");
    model.component("comp1").geom("geom1").feature("wp8").geom().feature("sq1").set("base", "center");
    model.component("comp1").geom("geom1").feature("wp8").geom().feature("sq1").set("size", 10);
    model.component("comp1").geom("geom1").create("ext5", "Extrude");
    model.component("comp1").geom("geom1").feature("ext5").set("workplane", "wp5");
    model.component("comp1").geom("geom1").feature("ext5").setIndex("distance", "10", 0);
    model.component("comp1").geom("geom1").feature("ext5").selection("input").set("wp5");
    model.component("comp1").geom("geom1").create("ext7", "Extrude");
    model.component("comp1").geom("geom1").feature("ext7").set("workplane", "wp7");
    model.component("comp1").geom("geom1").feature("ext7").setIndex("distance", "10", 0);
    model.component("comp1").geom("geom1").feature("ext7").selection("input").set("wp7");
    model.component("comp1").geom("geom1").create("ext8", "Extrude");
    model.component("comp1").geom("geom1").feature("ext8").setIndex("distance", "10", 0);
    model.component("comp1").geom("geom1").feature("ext8").selection("input").set("wp8");
    model.component("comp1").geom("geom1").create("ext6", "Extrude");
    model.component("comp1").geom("geom1").feature("ext6").set("workplane", "wp6");
    model.component("comp1").geom("geom1").feature("ext6").setIndex("distance", "10", 0);
    model.component("comp1").geom("geom1").feature("ext6").selection("input").set("wp6");
    model.component("comp1").geom("geom1").create("int1", "Intersection");
    model.component("comp1").geom("geom1").feature("int1").selection("input")
         .set("ext1", "ext2", "ext3", "ext4", "ext5", "ext6", "ext7", "ext8");
    model.component("comp1").geom("geom1").create("rt1", "RigidTransform");
    model.component("comp1").geom("geom1").feature("rt1").set("displ", new double[]{0, 0, -1.5});
    model.component("comp1").geom("geom1").feature("rt1").set("axis", new int[]{0, 0, 1});
    model.component("comp1").geom("geom1").feature("rt1").selection("input").set("int1");
    model.component("comp1").geom("geom1").create("blk1", "Block");
    model.component("comp1").geom("geom1").feature("blk1").set("pos", new int[]{0, 0, 0});
    model.component("comp1").geom("geom1").feature("blk1").set("base", "center");
    model.component("comp1").geom("geom1").feature("blk1").set("size", new int[]{20, 20, 20});
    model.component("comp1").geom("geom1").run();

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material().create("mat2", "Common");
    model.component("comp1").material("mat1").selection().set(2);
    model.component("comp1").material("mat2").selection().set(1);

    model.component("comp1").physics().create("es", "Electrostatics", "geom1");
    model.component("comp1").physics("es").create("pot1", "ElectricPotential", 2);
    model.component("comp1").physics("es").feature("pot1").selection().set(1, 2, 3, 4, 5, 14);

    model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");

    model.result().table("evl3").label("Evaluation 3D");
    model.result().table("evl3").comments("Interactive 3D values");

    model.component("comp1").view("view1").set("transparency", true);
    model.component("comp1").view("view2").axis().set("xmin", -5.714061737060547);
    model.component("comp1").view("view2").axis().set("xmax", 15.887458801269531);
    model.component("comp1").view("view2").axis().set("ymin", -10.442992210388184);
    model.component("comp1").view("view2").axis().set("ymax", 6.791007041931152);

    model.component("comp1").material("mat1").propertyGroup("def")
         .set("relpermittivity", new String[]{"-0.9135+4.9192i", "0", "0", "0", "-0.9135+4.9192i", "0", "0", "0", "-0.9135+4.9192i"});
    model.component("comp1").material("mat1").propertyGroup("def").set("relpermittivity_symmetry", "0");
    model.component("comp1").material("mat2").propertyGroup("def")
         .set("relpermittivity", new String[]{"1", "0", "0", "0", "1", "0", "0", "0", "1"});
    model.component("comp1").material("mat2").propertyGroup("def").set("relpermittivity_symmetry", "0");

    model.component("comp1").physics("es").feature("pot1").set("V0", "z[V/m]");

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
    model.result("pg1").create("arwl1", "ArrowLine");
    model.result("pg1").create("slc1", "Slice");
    model.result().export().create("data1", "Data");

    model.sol("sol1").attach("std1");
    model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("prefun", "amg");
    model.sol("sol1").runAll();

    model.result("pg1").label("Electric Potential (es)");
    model.result("pg1").set("frametype", "spatial");
    model.result("pg1").feature("mslc1").active(false);
    model.result("pg1").feature("mslc1").set("colortable", "RainbowLight");
    model.result("pg1").feature("mslc1").set("resolution", "normal");
    model.result("pg1").feature("arwl1").active(false);
    model.result("pg1").feature("arwl1").set("placement", "elements");
    model.result("pg1").feature("arwl1").set("scale", 0.4);
    model.result("pg1").feature("arwl1").set("scaleactive", true);
    model.result("pg1").feature("slc1").set("quickxnumber", 1);
    model.result("pg1").feature("slc1").set("resolution", "normal");
    model.result().export("data1").set("expr", new String[]{"V"});
    model.result().export("data1").set("unit", new String[]{"V"});
    model.result().export("data1").set("descr", new String[]{"Electric potential"});
    model.result().export("data1").set("filename", "/home/sjoao/octahedron/pot_freq2p7_rad8.0.dat");
    model.result().export("data1").set("location", "file");
    model.result().export("data1").set("coordfilename", "/home/sjoao/octahedron/positions_rad8.0.dat");
    model.result().export("data1").run();

    return model;
  }

  public static void main(String[] args) {
    run();
  }

}
