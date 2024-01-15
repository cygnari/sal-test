#include <vector>
#include <iostream>

#include "general_utils.hpp"
#include "structs.hpp"

void initialize_icosahedron(const RunConfig& run_information, std::vector<IcosPanel>& icos_panels, const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& z) {
  // sets up fast sum icosahedron
  icos_panels.resize(20);
  double phi = (1 + sqrt(5)) / 2;
  std::vector<double> point_1 = project_to_sphere(0, 1, phi, run_information.radius);
  std::vector<double> point_2 = project_to_sphere(0, -1, phi, run_information.radius);
  std::vector<double> point_3 = project_to_sphere(0, 1, -phi, run_information.radius);
  std::vector<double> point_4 = project_to_sphere(0, -1, -phi, run_information.radius);
  std::vector<double> point_5 = project_to_sphere(1, phi, 0, run_information.radius);
  std::vector<double> point_6 = project_to_sphere(1, -phi, 0, run_information.radius);
  std::vector<double> point_7 = project_to_sphere(-1, phi, 0, run_information.radius);
  std::vector<double> point_8 = project_to_sphere(-1, -phi, 0, run_information.radius);
  std::vector<double> point_9 = project_to_sphere(phi, 0, 1, run_information.radius);
  std::vector<double> point_10 = project_to_sphere(phi, 0, -1, run_information.radius);
  std::vector<double> point_11 = project_to_sphere(-phi, 0, 1, run_information.radius);
  std::vector<double> point_12 = project_to_sphere(-phi, 0, -1, run_information.radius);

  icos_panels[0].vertex_1 = point_1, icos_panels[0].vertex_2 = point_2, icos_panels[0].vertex_3 = point_9;
  icos_panels[1].vertex_1 = point_1, icos_panels[1].vertex_2 = point_2, icos_panels[1].vertex_3 = point_11;
  icos_panels[2].vertex_1 = point_1, icos_panels[2].vertex_2 = point_5, icos_panels[2].vertex_3 = point_7;
  icos_panels[3].vertex_1 = point_1, icos_panels[3].vertex_2 = point_5, icos_panels[3].vertex_3 = point_9;
  icos_panels[4].vertex_1 = point_1, icos_panels[4].vertex_2 = point_7, icos_panels[4].vertex_3 = point_11;
  icos_panels[5].vertex_1 = point_2, icos_panels[5].vertex_2 = point_6, icos_panels[5].vertex_3 = point_8;
  icos_panels[6].vertex_1 = point_2, icos_panels[6].vertex_2 = point_6, icos_panels[6].vertex_3 = point_9;
  icos_panels[7].vertex_1 = point_2, icos_panels[7].vertex_2 = point_8, icos_panels[7].vertex_3 = point_11;
  icos_panels[8].vertex_1 = point_3, icos_panels[8].vertex_2 = point_4, icos_panels[8].vertex_3 = point_10;
  icos_panels[9].vertex_1 = point_3, icos_panels[9].vertex_2 = point_4, icos_panels[9].vertex_3 = point_12;
  icos_panels[10].vertex_1 = point_3, icos_panels[10].vertex_2 = point_5, icos_panels[10].vertex_3 = point_7;
  icos_panels[11].vertex_1 = point_3, icos_panels[11].vertex_2 = point_5, icos_panels[11].vertex_3 = point_10;
  icos_panels[12].vertex_1 = point_3, icos_panels[12].vertex_2 = point_7, icos_panels[12].vertex_3 = point_12;
  icos_panels[13].vertex_1 = point_4, icos_panels[13].vertex_2 = point_6, icos_panels[13].vertex_3 = point_8;
  icos_panels[14].vertex_1 = point_4, icos_panels[14].vertex_2 = point_6, icos_panels[14].vertex_3 = point_10;
  icos_panels[15].vertex_1 = point_4, icos_panels[15].vertex_2 = point_8, icos_panels[15].vertex_3 = point_12;
  icos_panels[16].vertex_1 = point_5, icos_panels[16].vertex_2 = point_9, icos_panels[16].vertex_3 = point_10;
  icos_panels[17].vertex_1 = point_6, icos_panels[17].vertex_2 = point_9, icos_panels[17].vertex_3 = point_10;
  icos_panels[18].vertex_1 = point_7, icos_panels[18].vertex_2 = point_11, icos_panels[18].vertex_3 = point_12;
  icos_panels[19].vertex_1 = point_8, icos_panels[19].vertex_2 = point_11, icos_panels[19].vertex_3 = point_12;

  for (int i = 0; i < 20; i++) {
    icos_panels[i].id = i;
  }

  double xval, yval, zval;
  std::vector<double> points (9, 0), point (3, 0);
  bool point_found;
  int status;

  for (int i = 0; i < run_information.point_count; i++) {
    // assign the points to the base icosahedron
    xval = x[i];
    yval = y[i];
    zval = z[i];
    point_found = false;
    for (int j = 0; j < 20; j++) {
      point[0] = xval;
      point[1] = yval;
      point[2] = zval;
      points[0] = icos_panels[j].vertex_1[0];
      points[1] = icos_panels[j].vertex_1[1];
      points[2] = icos_panels[j].vertex_1[2];
      points[3] = icos_panels[j].vertex_2[0];
      points[4] = icos_panels[j].vertex_2[1];
      points[5] = icos_panels[j].vertex_2[2];
      points[6] = icos_panels[j].vertex_3[0];
      points[7] = icos_panels[j].vertex_3[1];
      points[8] = icos_panels[j].vertex_3[2];

      status = linear_solve(points, point, 3, 1, 1);
      if (status > 0) {
        throw std::runtime_error("Error with barycoordinate computationa in icosahedron initialize, line 77");
      }
      if ((point[0] >= 0) and (point[1] >= 0) and (point[2] >= 0)) {
        point_found = true;
        icos_panels[j].point_count++;
        icos_panels[j].points_inside.push_back(i);
      }

      if (point_found) {
        break;
      }
    }
    if (not point_found) {
      throw std::runtime_error("Not all points located");
    }
  }

  double x1, x2, x3, y1, y2, y3, z1, z2, z3, x12, x23, x31, y12, y23, y31, z12, z23, z31;
  std::vector<double> v12, v23, v31;
  int point_count;
  int start, end;
  double p_x, p_y, p_z;
  std::vector<int> points_to_assign;

  for (int i = 0; i < icos_panels.size(); i++) {
    if ((icos_panels[i].is_leaf) and (icos_panels[i].point_count > run_information.fast_sum_cluster_thresh)) {
      // refine, create 4 new sub-panels
      IcosPanel sub_panel1, sub_panel2, sub_panel3, sub_panel4;
      sub_panel1.parent_panel = &icos_panels[i];
      sub_panel2.parent_panel = &icos_panels[i];
      sub_panel3.parent_panel = &icos_panels[i];
      sub_panel4.parent_panel = &icos_panels[i];
      icos_panels[i].child_panel_1 = &sub_panel1;
      icos_panels[i].child_panel_2 = &sub_panel2;
      icos_panels[i].child_panel_3 = &sub_panel3;
      icos_panels[i].child_panel_4 = &sub_panel4;
      x1 = icos_panels[i].vertex_1[0], y1 = icos_panels[i].vertex_1[1], z1 = icos_panels[i].vertex_1[2];
      x2 = icos_panels[i].vertex_2[0], y2 = icos_panels[i].vertex_2[1], z2 = icos_panels[i].vertex_2[2];
      x3 = icos_panels[i].vertex_3[0], y3 = icos_panels[i].vertex_3[1], z3 = icos_panels[i].vertex_3[2];
      x12 = (x1 + x2) / 2, y12 = (y1 + y2) / 2, z12 = (z1 + z2) / 2;
      x23 = (x2 + x3) / 2, y23 = (y2 + y3) / 2, z23 = (z2 + z3) / 2;
      x31 = (x3 + x1) / 2, y31 = (y3 + y1) / 2, z31 = (z3 + z1) / 2;
      v12 = project_to_sphere(x12, y12, z12, run_information.radius);
      v23 = project_to_sphere(x23, y23, z23, run_information.radius);
      v31 = project_to_sphere(x31, y31, z31, run_information.radius);
      sub_panel1.vertex_1 = icos_panels[i].vertex_1, sub_panel1.vertex_2 = v31, sub_panel1.vertex_3 = v12;
      sub_panel2.vertex_1 = icos_panels[i].vertex_3, sub_panel2.vertex_2 = v23, sub_panel2.vertex_3 = v31;
      sub_panel3.vertex_1 = icos_panels[i].vertex_2, sub_panel3.vertex_2 = v12, sub_panel3.vertex_3 = v23;
      sub_panel4.vertex_1 = v12, sub_panel4.vertex_2 = v23, sub_panel4.vertex_3 = v31;
      sub_panel1.level = icos_panels[i].level + 1;
      sub_panel2.level = icos_panels[i].level + 1;
      sub_panel3.level = icos_panels[i].level + 1;
      sub_panel4.level = icos_panels[i].level + 1;
      icos_panels[i].is_leaf = false;
      start = icos_panels.size();
      sub_panel1.id = start;
      sub_panel2.id = start + 1;
      sub_panel3.id = start + 2;
      sub_panel4.id = start + 3;
      icos_panels.push_back(sub_panel1);
      icos_panels.push_back(sub_panel2);
      icos_panels.push_back(sub_panel3);
      icos_panels.push_back(sub_panel4);
      end = icos_panels.size();
      points_to_assign = icos_panels[i].points_inside;
      for (int j = 0; j < points_to_assign.size(); j++) {
        // loop over points in parent panel
        point_found = false;
        for (int k = start; k < end; k++) {
          // loop over child panels
          point[0] = x[points_to_assign[j]];
          point[1] = y[points_to_assign[j]];
          point[2] = z[points_to_assign[j]];
          points[0] = icos_panels[k].vertex_1[0];
          points[1] = icos_panels[k].vertex_1[1];
          points[2] = icos_panels[k].vertex_1[2];
          points[3] = icos_panels[k].vertex_2[0];
          points[4] = icos_panels[k].vertex_2[1];
          points[5] = icos_panels[k].vertex_2[2];
          points[6] = icos_panels[k].vertex_3[0];
          points[7] = icos_panels[k].vertex_3[1];
          points[8] = icos_panels[k].vertex_3[2];

          status = linear_solve(points, point, 3, 1, 1);
          if (status > 0) {
            throw std::runtime_error("Error with barycoordinate computation in icosahedron initialize, line 162");
          }

          if ((point[0] >= 0) and (point[1] >= 0) and (point[2] >= 0)) {
            point_found = true;
            icos_panels[k].point_count++;
            icos_panels[k].points_inside.push_back(points_to_assign[j]);
          }

          if (point_found) {
            break;
          }
        }
        if (not point_found) {
          throw std::runtime_error("Point not correctly found");
        }
      }
      point_count = icos_panels[start].point_count + icos_panels[start+1].point_count + icos_panels[start+2].point_count + icos_panels[start+3].point_count;
      if (point_count != icos_panels[i].point_count) {
        std::cout << "discrepancy: " << point_count << " vs " << icos_panels[i].point_count << std::endl;
        throw std::runtime_error("Not all points assigned");
      }
    }
  }

  double xc, yc, zc, d1, d2, d3, d;
  std::vector<double> vc;

  for (int i = 0; i < icos_panels.size(); i++) {
    x1 = icos_panels[i].vertex_1[0], y1 = icos_panels[i].vertex_1[1], z1 = icos_panels[i].vertex_1[2];
    x2 = icos_panels[i].vertex_2[0], y2 = icos_panels[i].vertex_2[1], z2 = icos_panels[i].vertex_2[2];
    x3 = icos_panels[i].vertex_3[0], y3 = icos_panels[i].vertex_3[1], z3 = icos_panels[i].vertex_3[2];
    xc = (x1 + x2 + x3) / 3.0;
    yc = (y1 + y2 + y3) / 3.0;
    zc = (z1 + z2 + z3) / 3.0;
    vc = project_to_sphere(xc, yc, zc, run_information.radius);
    xc = vc[0], yc = vc[1], zc = vc[2];
    icos_panels[i].center_p[0] = xc, icos_panels[i].center_p[1] = yc, icos_panels[i].center_p[2] = zc;
    d1 = pow(x1 - xc, 2) + pow(y1 - yc, 2) + pow(z1 - zc, 2);
    d2 = pow(x2 - xc, 2) + pow(y2 - yc, 2) + pow(z2 - zc, 2);
    d3 = pow(x3 - xc, 2) + pow(y3 - yc, 2) + pow(z3 - zc, 2);
    d = std::max(d1, std::max(d2, d3));
    icos_panels[i].radius = sqrt(d);
  }
}
