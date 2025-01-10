double get_polygon_area(Vector<Vector<double>> Polygon)
{
  double area = 0;
  unsigned n_vertex = Polygon.size();
  for (unsigned n = 0; n < n_vertex; n++)
  {
    double x0 = Polygon[n][0];
    double y0 = Polygon[n][1];
    double x1;
    double y1;
    if (n + 1 < n_vertex)
    {
      x1 = Polygon[n + 1][0];
      y1 = Polygon[n + 1][1];
    }
    else
    {
      x1 = Polygon[0][0];
      y1 = Polygon[0][1];
    }

    area += (x0 * y1 - y0 * x1) / 2;
  }
  return area;
}

void check_topology()
{
  double distance_tol = 1e-8;
  bool Check_for_pinchoff = true;
  bool Check_for_coalescence = true;

  bool Pinchoff_needed = false;

  for (unsigned i_bubble = 0;
       i_bubble < Problem_Parameter::Ordered_bubbles.size();
       i_bubble++)
  {
    double area =
      get_polygon_area(Problem_Parameter::Ordered_bubbles[i_bubble]);

    std::cout << "Bubble " << i_bubble << " has area " << area << std::endl;

    if (area < 0)
    {
      std::reverse(Problem_Parameter::Ordered_bubbles[i_bubble].begin(),
                   Problem_Parameter::Ordered_bubbles[i_bubble].end());
    }
  }


  if (Check_for_coalescence && Problem_Parameter::Ordered_bubbles.size() > 1)
  {
    std::cout << "Checking coalescence " << std::endl;
    for (unsigned i_detect = 0; i_detect < 5; i_detect++)
    {
      bool coalescence_needed = false;
      /// If there is more than one bubble, we can check the distance between
      /// them
      unsigned N_0 = Problem_Parameter::Ordered_bubbles.size();

      for (unsigned i_bubble = 0; i_bubble < N_0; i_bubble++)
      {
        for (unsigned j_bubble = 0; j_bubble < N_0; j_bubble++)
        {
          double approach_tol = 0.02;
          approach_tol = 10e-3;
          //	  approach_tol = 2e-4;
          //	  approach_tol = 0.01;
          //	  approach_tol = -1.0;
          if (i_bubble == j_bubble)
          {
            break;
          }
          /// Check the minimum approach distance between each pair of bubbles
          double min_distance = 100.0;
          double next_smallest = 100.0;
          unsigned i_vertex =
            Problem_Parameter::Ordered_bubbles[i_bubble].size();
          unsigned j_vertex =
            Problem_Parameter::Ordered_bubbles[j_bubble].size();
          for (unsigned i = 0; i < i_vertex; i++)
          {
            double x_i = Problem_Parameter::Ordered_bubbles[i_bubble][i][0];
            double y_i = Problem_Parameter::Ordered_bubbles[i_bubble][i][1];
            for (unsigned j = 0; j < j_vertex; j++)
            {
              double x_j = Problem_Parameter::Ordered_bubbles[j_bubble][j][0];
              double y_j = Problem_Parameter::Ordered_bubbles[j_bubble][j][1];
              if ((x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) <
                  min_distance)
              {
                next_smallest = min_distance;
                min_distance =
                  (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i);
              }
            }
          }
          min_distance = std::sqrt(min_distance);
          next_smallest = std::sqrt(next_smallest);
          std::cout << "Minimum distance between bubble " << i_bubble
                    << " and bubble " << j_bubble << " is " << min_distance
                    << std::endl;
          Problem_Parameter::Min_distance = min_distance;
          if (min_distance < approach_tol)
          {
            std::cout << "Setting approach tolerance to match minimum distance"
                      << std::endl;
            std::cout << "Measured min distance: " << min_distance
                      << " Second smallest: " << next_smallest
                      << " Current tolerance: " << approach_tol << std::endl;
            approach_tol = (min_distance + next_smallest) / 2.0;
            std::cout << "Set approach tolerance to " << approach_tol
                      << std::endl;
            std::cout << "Should mean we exclude only one node on each bubble "
                      << std::endl;

            coalescence_needed = true;
            Problem_Parameter::Topology_change_needed = true;
            std::cout << "Try to merge bubbles " << std::endl;

            /// Now need to merge bubble i and bubble j
            Vector<Vector<double>> Working_segment_vector;
            Vector<Vector<Vector<double>>> All_segments_vector;
            Vector<Vector<double>> Coalesced_i_vector;
            Vector<Vector<double>> Coalesced_j_vector;
            Vector<Vector<double>> Coalesced_vector;

            unsigned n_node =
              Problem_Parameter::Ordered_bubbles[i_bubble].size();
            std::cout << "Reject the following nodes on bubble i: "
                      << std::endl;

            for (unsigned i = 0; i < i_vertex; i++)
            {
              double x_i = Problem_Parameter::Ordered_bubbles[i_bubble][i][0];
              double y_i = Problem_Parameter::Ordered_bubbles[i_bubble][i][1];
              bool accept_node = true;
              /// Check that node i in bubble i_bubble is not too close to any
              /// node in bubble j_bubble
              for (unsigned j = 0; j < j_vertex; j++)
              {
                double x_j = Problem_Parameter::Ordered_bubbles[j_bubble][j][0];
                double y_j = Problem_Parameter::Ordered_bubbles[j_bubble][j][1];
                if ((x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) <
                    approach_tol * approach_tol)
                {
                  std::cout << x_i << " " << y_i << " "
                            << std::sqrt((x_j - x_i) * (x_j - x_i) +
                                         (y_j - y_i) * (y_j - y_i))
                            << std::endl;
                  accept_node = false;
                  break;
                }
              }
              if (accept_node == true)
              {
                /// If sufficiently far from all nodes on other bubbles, keep
                /// this node
                Working_segment_vector.push_back(
                  Problem_Parameter::Ordered_bubbles[i_bubble][i]);
              }
              else
              {
                /// If this is the first disallowed node, then we have a
                /// complete segment
                if (Working_segment_vector.size() > 0)
                {
                  All_segments_vector.push_back(Working_segment_vector);
                }
                Working_segment_vector.resize(0);
              }
            }
            /// Add the final segment if necessary
            if (Working_segment_vector.size() > 0)
            {
              All_segments_vector.push_back(Working_segment_vector);
            }
            Working_segment_vector.resize(0);

            std::cout << "Bubble " << i_bubble << " is in "
                      << All_segments_vector.size() << " segments "
                      << std::endl;

            /// Do we need to join two segments together?
            if (All_segments_vector.size() == 2)
            {
              /// Then join the second to the first
              std::cout << "Join first and last segments" << std::endl;
              Vector<Vector<double>> Joined_segment;
              n_node = All_segments_vector[1].size();
              for (unsigned i_node = 0; i_node < n_node; i_node++)
              {
                Joined_segment.push_back(All_segments_vector[1][i_node]);
              }
              n_node = All_segments_vector[0].size();
              for (unsigned i_node = 1; i_node < n_node; i_node++)
              {
                Joined_segment.push_back(All_segments_vector[0][i_node]);
              }
              All_segments_vector[0] = Joined_segment;
              All_segments_vector.erase(All_segments_vector.end());
              std::cout << "Successfully joined first and last segments"
                        << std::endl;
            }
            else if (All_segments_vector.size() > 2)
            {
              std::cout << "The bubble is trying to merge in too many segments "
                           "simultaneously"
                        << std::endl;
              std::cout << "Are you trying to enclose a drop!!!!?" << std::endl;
              std::cout << "Are you trying to enclose a drop!!!!?" << std::endl;
              for (unsigned i_segment = 0;
                   i_segment < All_segments_vector.size();
                   i_segment++)
              {
                std::cout << "Segment " << i_segment << " has size "
                          << All_segments_vector[i_segment].size() << std::endl;
                for (unsigned i_node = 0;
                     i_node < All_segments_vector[i_segment].size();
                     i_node++)
                {
                  std::cout << i_node << " "
                            << All_segments_vector[i_segment][i_node][0] << " "
                            << All_segments_vector[i_segment][i_node][1]
                            << std::endl;
                }
              }
            }

            Coalesced_i_vector = All_segments_vector[0];

            All_segments_vector.resize(0);
            n_node = Problem_Parameter::Ordered_bubbles[j_bubble].size();

            std::cout << "Reject the following nodes on bubble j " << std::endl;
            for (unsigned j = 0; j < j_vertex; j++)
            {
              double x_j = Problem_Parameter::Ordered_bubbles[j_bubble][j][0];
              double y_j = Problem_Parameter::Ordered_bubbles[j_bubble][j][1];
              bool accept_node = true;

              for (unsigned i = 0; i < i_vertex; i++)
              {
                double x_i = Problem_Parameter::Ordered_bubbles[i_bubble][i][0];
                double y_i = Problem_Parameter::Ordered_bubbles[i_bubble][i][1];
                if ((x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i) <
                    approach_tol * approach_tol)
                {
                  std::cout << x_j << " " << y_j << " "
                            << std::sqrt((x_j - x_i) * (x_j - x_i) +
                                         (y_j - y_i) * (y_j - y_i))
                            << std::endl;

                  accept_node = false;
                  break;
                }
              }
              if (accept_node == true)
              {
                Working_segment_vector.push_back(
                  Problem_Parameter::Ordered_bubbles[j_bubble][j]);
              }
              else
              {
                if (Working_segment_vector.size() > 0)
                {
                  All_segments_vector.push_back(Working_segment_vector);
                }
                Working_segment_vector.resize(0);
              }
            }
            if (Working_segment_vector.size() > 0)
            {
              All_segments_vector.push_back(Working_segment_vector);
            }
            Working_segment_vector.resize(0);

            std::cout << "Bubble " << j_bubble << " is in "
                      << All_segments_vector.size() << " segments "
                      << std::endl;

            if (All_segments_vector.size() == 2)
            {
              /// Then join the second to the first
              std::cout << "Join first and last segments" << std::endl;
              Vector<Vector<double>> Joined_segment;
              n_node = All_segments_vector[1].size();
              for (unsigned i_node = 0; i_node < n_node; i_node++)
              {
                Joined_segment.push_back(All_segments_vector[1][i_node]);
              }
              n_node = All_segments_vector[0].size();
              for (unsigned i_node = 1; i_node < n_node; i_node++)
              {
                Joined_segment.push_back(All_segments_vector[0][i_node]);
              }
              All_segments_vector[0] = Joined_segment;
              All_segments_vector.erase(All_segments_vector.end());
              std::cout << "Successfully joined first and last segments"
                        << std::endl;
            }
            else if (All_segments_vector.size() > 2)
            {
              std::cout << "The bubble is trying to merge in too many segments "
                           "simultaneously"
                        << std::endl;
              std::cout << "Are you trying to enclose a drop!!!!?" << std::endl;
              for (unsigned i_segment = 0;
                   i_segment < All_segments_vector.size();
                   i_segment++)
              {
                std::cout << "Segment " << i_segment << " has size "
                          << All_segments_vector[i_segment].size() << std::endl;
                for (unsigned i_node = 0;
                     i_node < All_segments_vector[i_segment].size();
                     i_node++)
                {
                  std::cout << i_node << " "
                            << All_segments_vector[i_segment][i_node][0] << " "
                            << All_segments_vector[i_segment][i_node][1]
                            << std::endl;
                }
              }
            }
            Coalesced_j_vector = All_segments_vector[0];

            /// Finally, join parts i and j, in any order
            n_node = Coalesced_i_vector.size();
            for (unsigned i_node = 0; i_node < n_node; i_node++)
            {
              Coalesced_vector.push_back(Coalesced_i_vector[i_node]);
            }
            n_node = Coalesced_j_vector.size();
            for (unsigned j_node = 0; j_node < n_node; j_node++)
            {
              Coalesced_vector.push_back(Coalesced_j_vector[j_node]);
            }
            /// Close the loop
            Coalesced_vector.push_back(Coalesced_i_vector[0]);

            /// Replace the first bubble with the new coalesced version
            Problem_Parameter::Ordered_bubbles[i_bubble] = Coalesced_vector;
            /// And remove the second bubble
            Problem_Parameter::Ordered_bubbles.erase(
              Problem_Parameter::Ordered_bubbles.begin() + j_bubble);
          } // If min distance
          if (coalescence_needed == true)
          {
            break;
          }
        } // Test bubble j
        if (coalescence_needed == true)
        {
          break;
        }
      } // Test bubble i
      if (coalescence_needed == false)
      {
        break;
      }
      else
      {
        std::cout << "Check neighbouring points " << std::endl;
        /// Remove any points which are too close to their neighbours
        double neighbouring_distance_tol = 1e-4;
        for (unsigned bubble = 0;
             bubble < Problem_Parameter::Ordered_bubbles.size();
             bubble++)
        {
          unsigned n_node = Problem_Parameter::Ordered_bubbles[bubble].size();
          double x0 = Problem_Parameter::Ordered_bubbles[bubble][0][0];
          double y0 = Problem_Parameter::Ordered_bubbles[bubble][0][1];
          Vector<Vector<double>> Accepted_vertices;
          Accepted_vertices.push_back(
            Problem_Parameter::Ordered_bubbles[bubble][0]);
          for (unsigned i_vertex = 1; i_vertex < n_node; i_vertex++)
          {
            double x1 = Problem_Parameter::Ordered_bubbles[bubble][i_vertex][0];
            double y1 = Problem_Parameter::Ordered_bubbles[bubble][i_vertex][1];
            if ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) >
                neighbouring_distance_tol)
            {
              Accepted_vertices.push_back(
                Problem_Parameter::Ordered_bubbles[bubble][i_vertex]);
              x0 = x1;
              y0 = y1;
            }
            else
            {
              std::cout << "Removed node as too close" << std::endl;
              std::cout << x0 << " " << y0 << " " << x1 << " " << y1 << " "
                        << (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)
                        << std::endl;
            }
          }

          /// Check that first and last vertices are the same
          double x1 = Problem_Parameter::Ordered_bubbles[bubble][0][0];
          double y1 = Problem_Parameter::Ordered_bubbles[bubble][0][1];
          x0 = Problem_Parameter::Ordered_bubbles[bubble][n_node - 1][0];
          y0 = Problem_Parameter::Ordered_bubbles[bubble][n_node - 1][1];
          if ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) > 1e-12)
          {
            Accepted_vertices.push_back(
              Problem_Parameter::Ordered_bubbles[bubble][0]);
            std::cout << "Add another copy of the first node " << std::endl;
          }

          std::cout << "Bubble " << bubble << " had size "
                    << Problem_Parameter::Ordered_bubbles[bubble].size()
                    << " and now has size " << Accepted_vertices.size()
                    << std::endl;
          Problem_Parameter::Ordered_bubbles[bubble] = Accepted_vertices;
        }
      }
    }
  }

  if (Check_for_pinchoff)
  {
    /// Each bubble only needs to be checked once
    for (unsigned bubble = 0;
         bubble < Problem_Parameter::Ordered_bubbles.size();
         bubble++)
    {
      unsigned n_node = Problem_Parameter::Ordered_bubbles[bubble].size();
      Vector<unsigned> Allowable_node(n_node, 1);

      Vector<Vector<Vector<double>>> All_segments_vector;
      Vector<Vector<double>> Working_segment_vector;
      unsigned Working_segment_label = 0;

      double initial_area =
        get_polygon_area(Problem_Parameter::Ordered_bubbles[bubble]);

      Pinchoff_needed = false;

      std::cout
        << "Check for segments crossing each other, and discard if necessary "
        << std::endl;

      unsigned Collision_number = 2;


      for (unsigned i = 0; i < n_node - 1; i++)
      {
        double x_a = Problem_Parameter::Ordered_bubbles[bubble][i][0];
        double y_a = Problem_Parameter::Ordered_bubbles[bubble][i][1];
        double x_b = Problem_Parameter::Ordered_bubbles[bubble][i + 1][0];
        double y_b = Problem_Parameter::Ordered_bubbles[bubble][i + 1][1];

        double na_x = y_b - y_a;
        double na_y = x_a - x_b;

        for (unsigned j = 0; j < n_node - 1; j++)
        {
          /// We are not worried about neighbouring segments
          if ((i != j && j != i + 1 && i != j + 1) &&
              (i != 0 || j != n_node - 2))
          {
            double x_0 = Problem_Parameter::Ordered_bubbles[bubble][j][0];
            double y_0 = Problem_Parameter::Ordered_bubbles[bubble][j][1];
            double x_1 = Problem_Parameter::Ordered_bubbles[bubble][j + 1][0];
            double y_1 = Problem_Parameter::Ordered_bubbles[bubble][j + 1][1];

            double n0_x = y_1 - y_0;
            double n0_y = x_0 - x_1;

            double mu = (na_x * (x_b - x_1) + na_y * (y_b - y_1)) /
                        (na_x * (x_0 - x_1) + na_y * (y_0 - y_1));
            double lambda = (n0_x * (x_1 - x_b) + n0_y * (y_1 - y_b)) /
                            (n0_x * (x_a - x_b) + n0_y * (y_a - y_b));
            // 	      std::cout << mu << " " << lambda << std::endl;

            if (mu > 0 && mu < 1.0 && lambda > 0 && lambda < 1.0)
            {
              Pinchoff_needed = true;
              std::cout << "Current status: " << Allowable_node[i] << " "
                        << Allowable_node[i + 1] << " " << Allowable_node[j]
                        << " " << Allowable_node[j + 1] << std::endl;
              if (Allowable_node[i] == 1)
              {
                Allowable_node[i] = Collision_number;
              }
              if (Allowable_node[i + 1] == 1)
              {
                Allowable_node[i + 1] = Collision_number;
              }
              if (Allowable_node[j] == 1)
              {
                Allowable_node[j] = Collision_number;
              }
              if (Allowable_node[j + 1] == 1)
              {
                Allowable_node[j + 1] = Collision_number;
              }
              std::cout << "Update to: " << Allowable_node[i] << " "
                        << Allowable_node[i + 1] << " " << Allowable_node[j]
                        << " " << Allowable_node[j + 1] << std::endl;
              Collision_number++;
              std::cout << " ij: " << i << " " << j << std::endl;
              /// The segment between nodes i and i+1 has an intersection;
              /// eliminate both nodes?
              break;
            }
          }
        }
      }
      if (Collision_number == 2)
      {
        std::cout << "No intersections" << std::endl;
        /// There are no intersections!
        continue;
        //                break;
      }
      else
      {
        unsigned Start_label = 2;
        for (unsigned i = 0; i < n_node - 1; i++)
        {
          std::cout << Problem_Parameter::Ordered_bubbles[bubble][i][0] << " "
                    << Problem_Parameter::Ordered_bubbles[bubble][i][1] << " "
                    << Allowable_node[i] << std::endl;

          if (Allowable_node[i] == 1)
          {
            Working_segment_vector.push_back(
              Problem_Parameter::Ordered_bubbles[bubble][i]);
          }
          else if (Allowable_node[i - 1] == 1)
          {
            unsigned End_label = Allowable_node[i];

            if (End_label == Start_label)
            {
              std::cout << "Retaining segment " << std::endl;
              /// We need to tie off a segment
              All_segments_vector.push_back(Working_segment_vector);
            }
            else
            {
              std::cout << "Discarding segment" << std::endl;
              // The segment connects two different intersections, and cannot
              // close.
            }
            Working_segment_vector.resize(0);
            Start_label = End_label;
          }
        }
      }

      /// We run out of points for the last segment, so may need to add
      /// manually.
      if (Working_segment_vector.size() > 0)
      {
        All_segments_vector.push_back(Working_segment_vector);
        Working_segment_label++;
      }
      // 	End_collision_number.push_back(0);
      /// If last node of last segment is the same as the first node of the
      /// first segment, then join (without repeating the node)

      if (Pinchoff_needed == true)
      {
        Problem_Parameter::Topology_change_needed = true;

        double x0 = All_segments_vector[0][0][0];
        double y0 = All_segments_vector[0][0][1];

        unsigned n_segments = All_segments_vector.size();
        unsigned n_node = All_segments_vector[n_segments - 1].size();
        double x1 = All_segments_vector[n_segments - 1][n_node - 1][0];
        double y1 = All_segments_vector[n_segments - 1][n_node - 1][1];

        std::cout << x0 << " " << y0 << " compare to " << x1 << " " << y1
                  << std::endl;

        if ((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1) < 0.1)
        {
          std::cout << "Join first and last segments" << std::endl;
          Vector<Vector<double>> Joined_segment;
          for (unsigned i_node = 0; i_node < n_node; i_node++)
          {
            Joined_segment.push_back(
              All_segments_vector[n_segments - 1][i_node]);
          }
          n_node = All_segments_vector[0].size();
          for (unsigned i_node = 1; i_node < n_node; i_node++)
          {
            Joined_segment.push_back(All_segments_vector[0][i_node]);
          }
          All_segments_vector[0] = Joined_segment;
          All_segments_vector.erase(All_segments_vector.end());
          std::cout << "Successfully joined first and last segments"
                    << std::endl;
        }


        std::cout << "Bubble has been divided into "
                  << All_segments_vector.size() << " segments" << std::endl;

        /// Remove any under-size segments
        n_segments = All_segments_vector.size();
        for (unsigned j = 0; j < n_segments; j++)
        {
          double area = get_polygon_area(All_segments_vector[j]);
          std::cout << "Segment " << j << " has area " << area << std::endl;
        }


        Vector<Vector<Vector<double>>> Closed_segments;
        Vector<Vector<double>> Misc_segments;
        Vector<unsigned> Segments_assigned(All_segments_vector.size(), 0);
        for (unsigned s = 0; s < All_segments_vector.size(); s++)
        {
          double segment_area = get_polygon_area(All_segments_vector[s]);
          std::cout << "Area for segment " << s << " is " << segment_area
                    << std::endl;

          //                     unsigned n_vertex =
          //                     All_segments_vector[s].size(); double x0 =
          //                     All_segments_vector[s][0][0]; double y0 =
          //                     All_segments_vector[s][0][1]; double x1 =
          //                     All_segments_vector[s][n_vertex-1][0]; double
          //                     y1 = All_segments_vector[s][n_vertex-1][1];
          //         //            double end_distance =
          //         std::sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));

          /// Keep only segments in the same sense as the original
          if (initial_area * segment_area > 0)
          {
            All_segments_vector[s].push_back(All_segments_vector[s][0]);
            Closed_segments.push_back(All_segments_vector[s]);
            Segments_assigned[s] = 1;
          }
        } // End loop of segments
        // 	  if (Misc_segments.size()>0)
        // 	  {
        // 	    Misc_segments.push_back(Misc_segments[0]);
        // 	    Closed_segments.push_back(Misc_segments);
        // 	  }


        std::cout << "Check neighbouring points " << std::endl;
        /// Remove any points which are too close to their neighbours
        distance_tol = 1e-3;
        for (unsigned segment = 0; segment < Closed_segments.size(); segment++)
        {
          unsigned n_node = Closed_segments[segment].size();
          std::cout << "Segment " << segment << " has " << n_node << " nodes"
                    << std::endl;
          double x0 = Closed_segments[segment][0][0];
          double y0 = Closed_segments[segment][0][1];
          Vector<Vector<double>> Accepted_vertices;
          Accepted_vertices.push_back(Closed_segments[segment][0]);
          for (unsigned i_vertex = 1; i_vertex < n_node; i_vertex++)
          {
            double x1 = Closed_segments[segment][i_vertex][0];
            double y1 = Closed_segments[segment][i_vertex][1];
            if ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) >
                distance_tol * distance_tol)
            {
              Accepted_vertices.push_back(Closed_segments[segment][i_vertex]);
            }
            x0 = x1;
            y0 = y1;
          }

          /// Check that first and last vertices are the same, and fix if needed
          double x1 = Closed_segments[segment][0][0];
          double y1 = Closed_segments[segment][0][1];
          if ((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) >
              distance_tol * distance_tol)
          {
            std::cout << "Add an extra copy of the first node" << std::endl;
            Accepted_vertices.push_back(Closed_segments[segment][0]);
          }

          std::cout << "Segment " << segment << " had size "
                    << Closed_segments[segment].size() << " and now has size "
                    << Accepted_vertices.size() << std::endl;
          Closed_segments[segment] = Accepted_vertices;
        }

        std::cout << "There are " << Closed_segments.size()
                  << " closed segments" << std::endl;
        std::cout << "Replace original bubble" << std::endl;
        Problem_Parameter::Ordered_bubbles[bubble] = Closed_segments[0];

        for (unsigned segment = 1; segment < Closed_segments.size(); segment++)
        {
          std::cout << "And append extra bubbles " << std::endl;
          Problem_Parameter::Ordered_bubbles.push_back(
            Closed_segments[segment]);
        }

        /// Check that each bubble has positive area.
        for (unsigned i_bubble = 0;
             i_bubble < Problem_Parameter::Ordered_bubbles.size();
             i_bubble++)
        {
          double area =
            get_polygon_area(Problem_Parameter::Ordered_bubbles[i_bubble]);
          if (initial_area * area < 0)
          {
            std::cout << "Need to remove this bubble " << std::endl;
            Problem_Parameter::Ordered_bubbles.erase(
              Problem_Parameter::Ordered_bubbles.begin() + i_bubble);
          }
        }
      } //// Was a topology change needed?
      std::cout << "Finished checking topology in bubble " << bubble
                << std::endl;
    } /// Loop over all the bubbles that originally existed.... TBH, this will
      /// go wrong if there's more than 1.
  }

  if (Problem_Parameter::Topology_change_needed == true)
  {
    Problem_Parameter::All_the_bubbles = Problem_Parameter::Ordered_bubbles;
    std::cout << "A topology change is required " << std::endl;
    Problem_Parameter::N_Bubble = Problem_Parameter::All_the_bubbles.size();
    Problem_Parameter::Reload_from_vector = true;
    std::cout << "Will try to construct " << Problem_Parameter::N_Bubble
              << " bubbles from vector " << std::endl;
  }

  /// Check for close approach


  std::cout << "Have finished ordering" << std::endl;
}
