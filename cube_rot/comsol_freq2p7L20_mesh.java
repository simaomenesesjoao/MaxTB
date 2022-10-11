/*
 * cube_rot.java
 */

import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on Oct 10 2022, 12:10 by COMSOL 5.5.0.359. */
public class comsol_freq2p7L20_mesh {

  public static Model run() {
    Model model = ModelUtil.create("Model");

    model.modelPath("/data/users/sjoao/cube_rot");

    model.label("cube_rot.mph");

    model.component().create("comp1", true);

    model.component("comp1").geom().create("geom1", 3);

    model.result().table().create("evl3", "Table");

    model.component("comp1").mesh().create("mesh1");

    model.component("comp1").geom("geom1").lengthUnit("nm");
    model.component("comp1").geom("geom1").create("blk1", "Block");
    model.component("comp1").geom("geom1").feature("blk1").set("pos", new int[]{0, 0, 0});
    model.component("comp1").geom("geom1").feature("blk1").set("base", "center");
    model.component("comp1").geom("geom1").feature("blk1").set("size", new int[]{20, 20, 20});
    model.component("comp1").geom("geom1").create("blk2", "Block");
    model.component("comp1").geom("geom1").feature("blk2").set("base", "center");
    model.component("comp1").geom("geom1").feature("blk2").set("size", new int[]{3, 3, 3});
    model.component("comp1").geom("geom1").create("rot1", "Rotate");
    model.component("comp1").geom("geom1").feature("rot1").set("axis", new int[]{1, 0, 0});
    model.component("comp1").geom("geom1").feature("rot1").setIndex("rot", "45", 0);
    model.component("comp1").geom("geom1").feature("rot1").selection("input").set("blk2");
    model.component("comp1").geom("geom1").create("sca1", "Scale");
    model.component("comp1").geom("geom1").feature("sca1").set("factor", 0.7071067);
    model.component("comp1").geom("geom1").feature("sca1").selection("input").set("rot1");
    model.component("comp1").geom("geom1").run();

    model.component("comp1").material().create("mat1", "Common");
    model.component("comp1").material().create("mat2", "Common");
    model.component("comp1").material("mat1").selection().set(2);
    model.component("comp1").material("mat2").selection().set(1);

    model.component("comp1").physics().create("es", "Electrostatics", "geom1");
    model.component("comp1").physics("es").create("pot1", "ElectricPotential", 2);
    model.component("comp1").physics("es").feature("pot1").selection().set(1, 2, 3, 4, 5, 12);

    model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");

    model.result().table("evl3").label("Evaluation 3D");
    model.result().table("evl3").comments("Interactive 3D values");

    model.component("comp1").view("view1").set("transparency", true);

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
    model.result().export("data1").set("filename", "/home/sjoao/data/cube_rot/pot_freq2p7L20_mesh.dat");
    model.result().export("data1").set("location", "file");
    model.result().export("data1").set("coordfilename", "/home/sjoao/data/cube_rot/positions_mesh.dat");
    model.result().export("data1").run();

    return model;
  }

  public static void main(String[] args) {
    run();
  }

}
