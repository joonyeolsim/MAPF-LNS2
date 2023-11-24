#include "SIPP.h"

void SIPP::updatePath(const LLNode* goal, Path& path) {
  num_collisions = goal->num_of_conflicts;
  path.resize(goal->state.timestep + 1);
  // num_of_conflicts = goal->num_of_conflicts;

  const auto* curr = goal;
  while (curr->parent != nullptr)  // non-root node
  {
    const auto* prev = curr->parent;
    int t = prev->state.timestep + 1;
    while (t < curr->state.timestep) {
      path[t].location = prev->state.location;  // wait at prev location
      path[t].timestep = t;
      path[t].orientation = prev->state.orientation;
      t++;
    }
    path[curr->state.timestep].location = curr->state.location;  // move to curr location
    path[curr->state.timestep].timestep = curr->state.timestep;
    path[curr->state.timestep].orientation = curr->state.orientation;
    curr = prev;
  }
  assert(curr->state.timestep == 0);
  path[0].location = curr->state.location;
  path[0].timestep = curr->state.timestep;
  path[0].orientation = curr->state.orientation;
}

// find path by A*
// Returns a path that minimizes the collisions with the paths in the path table, breaking ties by
// the length
Path SIPP::findPath(const ConstraintTable& constraint_table, vector<Path>& paths) {
  reset();
  ReservationTable reservation_table(constraint_table, goal_location, instance.window);
  Path path;
  Interval interval = reservation_table.get_first_safe_interval(start_state.location);
  if (get<0>(interval) > 0) return path;
  auto holding_time = constraint_table.getHoldingTime(goal_location, constraint_table.length_min);
  auto last_target_collision_time = constraint_table.getLastCollisionTimestep(goal_location);
  // generate start and add it to the OPEN & FOCAL list
  auto h = max(max(my_heuristic[start_state.location], holding_time), last_target_collision_time + 1);
  auto start = new SIPPNode(State(start_state.location, 0, start_state.orientation), 0, h, nullptr, get<1>(interval),
                            get<1>(interval), get<2>(interval), get<2>(interval));
 pushNodeToFocal(start);
  vector<tuple<int, int, int, int>> disfavor_table; // from, to, from_time, to_time
  vector<tuple<int, int, int, int, int>> prefer_table; // from, from_ori, to_ori, from_time, to_time

  while (!focal_list.empty()) {
    auto* curr = focal_list.top();
    focal_list.pop();
    curr->in_openlist = false;
    num_expanded++;
    assert(curr->state.location >= 0);
    // check if the popped node is a goal
    if (curr->is_goal) {
      updatePath(curr, path);
      break;
    } else if (curr->state.location == goal_location &&  // arrive at the goal location
               !curr->wait_at_goal &&                    // not wait at the goal location
               curr->state.timestep >=
                   holding_time)  // the agent can hold the goal location afterward
    {
      int future_collisions =
          constraint_table.getFutureNumOfCollisions(curr->state.location, curr->state.timestep);
      if (future_collisions == 0) {
        updatePath(curr, path);
        break;
      }
      // generate a goal node
      auto goal = new SIPPNode(*curr);
      goal->is_goal = true;
      goal->h_val = 0;
      goal->num_of_conflicts += future_collisions;
      // try to retrieve it from the hash table
      if (dominanceCheck(goal))
        pushNodeToFocal(goal);
      else
        delete goal;
    }

    for (auto& next_state : instance.getNeighbors(curr->state))  // move to neighboring locations
    {
      for (auto& i : reservation_table.get_safe_intervals(curr->state.location, next_state.location,
                                                          curr->state.timestep + 1,
                                                          curr->high_expansion + 1)) {
        int next_high_generation, next_timestep, next_high_expansion;
        int next_v_collision, next_e_collision;
        tie(next_high_generation, next_timestep, next_high_expansion, next_v_collision,
            next_e_collision) = i;

        if (next_v_collision > 0 || next_e_collision > 0) {
                  // Deadlock Detection
        auto &table = reservation_table.constraint_table.path_table_for_CAT->table;
        auto to = next_state.location;
        auto from = curr->state.location;
        auto to_time = next_timestep;
        auto from_time = to_time - 1;

        // scoring 1
        if (from != to && to_time - from_time == 1 && reservation_table.constraint_table.path_table_for_CAT != nullptr && !table.empty()) { // next state is moving forward
          if (table[to].size() > from_time && table[from].size() > to_time) {
            for (auto a1 : table[to][from_time]) {
              for (auto a2 : table[from][to_time]) {
                if (a1 == a2) { // the agent in front move forward
                  int degree = abs(paths[a1][from_time].orientation - curr->state.orientation);
                  if (degree == 2) { // the agents are facing each other
                    disfavor_table.emplace_back(from, to, from_time, to_time);
                    if (to_time + 1 < instance.window)
                      disfavor_table.emplace_back(from, to, from_time + 1, to_time + 1);
                    // cout << "Score1: " << "disfavor: " << from << " " << to << " " << from_time << " " << to_time << " " << degree << endl;
                    // cout << "Score1: " << "disfavor: " << from << " " << to << " " << from_time + 1 << " " << to_time + 1 << " " << degree << endl;
                  }
                }
              }
            }
          }

          // scoring 2
          if (table[to].size() > from_time + 1 && table[from].size() > to_time + 1) {
            for (auto a1 : table[to][from_time + 1]) {
              for (auto a2 : table[from][to_time + 1]) {
                if (a1 == a2) { // the agent in front wait and move forward
                  int degree = abs(paths[a1][from_time].orientation - curr->state.orientation);
                  if (degree == 1 || degree == 3 && paths[a1][from_time].location == to) { // the agents will face each other next step
                    disfavor_table.emplace_back(from, to, from_time, to_time);
                    // cout << "Score2: " << "disfavor: " << from << " " << to << " " << from_time << " " << to_time << " " << degree << endl;
                  }
                  if (degree == 2 && paths[a1][from_time].location == to) { // the agents will face each other next step
                    disfavor_table.emplace_back(from, to, from_time, to_time);
                    if (to_time + 1 < instance.window)
                      disfavor_table.emplace_back(from, to, from_time + 1, to_time + 1);
                    // cout << "Score2: " << "disfavor: " << from << " " << to << " " << from_time << " " << to_time << " " << degree << endl;
                    // cout << "Score2: " << "disfavor: " << from << " " << to << " " << from_time + 1 << " " << to_time + 1 << " " << degree << endl;
                  }
                }
              }
            }
          }

          // deadlock resolving 1
          if (table[to].size() > from_time + 2 && table[from].size() > to_time + 2) {
            for (auto a1 : table[to][from_time + 2]) {
              for (auto a2 : table[from][to_time + 2]) {
                if (a1 == a2 &&
                  paths[a1][from_time].location == to &&
                  paths[a1][from_time + 1].location == to) { // the agent in front wait and wait and move forward
                  assert((from_time == 0 && to_time == 1) || (from_time == 1 && to_time == 2));
                  int prefer_orientation;
                  // left rotation
                  if (curr->state.orientation == 0)
                    prefer_orientation = 3;
                  else
                    prefer_orientation = curr->state.orientation - 1;
                  prefer_table.emplace_back(curr->state.location, curr->state.orientation, prefer_orientation, from_time, to_time);
                  // cout << "prefer: " << curr->state.location << " " << curr->state.orientation << " " << prefer_orientation << " " << from_time << " " << to_time << endl;

                  // right rotation
                  if (curr->state.orientation == 3)
                    prefer_orientation = 0;
                  else
                    prefer_orientation = curr->state.orientation + 1;
                  prefer_table.emplace_back(curr->state.location, curr->state.orientation, prefer_orientation, from_time, to_time);
                  // cout << "prefer: " << curr->state.location << " " << curr->state.orientation << " " << prefer_orientation << " " << from_time << " " << to_time << endl;
                }
              }
            }
          }

          // deadlock resolving 1
          if (table[to].size() > from_time + 3 && table[from].size() > to_time + 3) {
            for (auto a1 : table[to][from_time + 3]) {
              for (auto a2 : table[from][to_time + 3]) {
                if (a1 == a2 &&
                  paths[a1][from_time].location == to &&
                  paths[a1][from_time + 1].location == to &&
                  paths[a1][from_time + 2].location == to) { // the agent in front wait and wait and move forward
                  assert(from_time == 0 && to_time == 1);
                  int prefer_orientation;
                  // left rotation
                  if (curr->state.orientation == 0)
                    prefer_orientation = 3;
                  else
                    prefer_orientation = curr->state.orientation - 1;
                  prefer_table.emplace_back(curr->state.location, curr->state.orientation, prefer_orientation, from_time, to_time);
                  // cout << "prefer: " << curr->state.location << " " << curr->state.orientation << " " << prefer_orientation << " " << from_time << " " << to_time << endl;

                  // right rotation
                  if (curr->state.orientation == 3)
                    prefer_orientation = 0;
                  else
                    prefer_orientation = curr->state.orientation + 1;
                  prefer_table.emplace_back(curr->state.location, curr->state.orientation, prefer_orientation, from_time, to_time);
                  // cout << "prefer: " << curr->state.location << " " << curr->state.orientation << " " << prefer_orientation << " " << from_time << " " << to_time << endl;
                }
              }
            }
          }
        }
        }

        // deadlock resolving 2
        if (!disfavor_table.empty()) {
          for (auto &disfavor : disfavor_table) {
            int from, to, from_time, to_time;
            tie(from, to, from_time, to_time) = disfavor;
            if (
              curr->state.location == from &&
              next_state.location == to &&
              next_timestep - 1 == from_time &&
              next_timestep == to_time
              ) {
              next_e_collision = INT_MAX / 2;
              // remove the score from the scoring table
              disfavor_table.erase(remove(disfavor_table.begin(), disfavor_table.end(), disfavor), disfavor_table.end());
            }
          }
        }

        // deadlock resolving 3
        if (!prefer_table.empty()) {
          for (auto &prefer : prefer_table) {
            int from, from_ori, to_ori, from_time, to_time;
            tie(from, from_ori, to_ori, from_time, to_time) = prefer;
            if (
              curr->state.location == from &&
              curr->state.orientation == from_ori &&
              next_state.orientation == to_ori &&
              next_timestep - 1 == from_time &&
              next_timestep == to_time
              ) {
              next_e_collision = -2;
              // remove the score from the scoring table
              prefer_table.erase(remove(prefer_table.begin(), prefer_table.end(), prefer), prefer_table.end());
              }
          }
        }

        if (next_timestep + my_heuristic[next_state.location] > constraint_table.length_max) break;
        auto next_collisions = curr->num_of_conflicts + next_v_collision + next_e_collision;
        auto next_h_val =
            max(my_heuristic[next_state.location],
                (next_collisions > 0 ? holding_time : curr->getFVal()) - next_timestep);
        // path max
        // generate (maybe temporary) node
        auto next = new SIPPNode(State(next_state.location, next_timestep, next_state.orientation),
                                 next_timestep, next_h_val, curr, next_high_generation,
                                 next_high_expansion, next_v_collision, next_collisions);
        // try to retrieve it from the hash table
        if (dominanceCheck(next))
          pushNodeToFocal(next);
        else
          delete next;
      }
    }  // end for loop that generates successors
    // wait at the current location
    if (curr->high_expansion == curr->high_generation and
        reservation_table.find_safe_interval(interval, curr->state.location,
                                             curr->high_expansion) and
        get<0>(interval) + curr->h_val <= reservation_table.constraint_table.length_max) {
      auto next_timestep = get<0>(interval);
      auto next_h_val =
          max(my_heuristic[curr->state.location],
              (get<2>(interval) ? holding_time : curr->getFVal()) - next_timestep);  // path max
      auto next_collisions =
          curr->num_of_conflicts +
          // (int)curr->collision_v * max(next_timestep - curr->state.timestep - 1, 0) +
          (int)get<2>(interval);
      auto next = new SIPPNode(State(curr->state.location, next_timestep, curr->state.orientation),
                               next_timestep, next_h_val, curr, get<1>(interval), get<1>(interval),
                               get<2>(interval), next_collisions);
      next->wait_at_goal = (curr->state.location == goal_location);
      if (dominanceCheck(next))
        pushNodeToFocal(next);
      else
        delete next;
    }
  }  // end while loop

  releaseNodes();
  return path;
}

// TODO:: currently this is implemented in SIPP inefficiently
int SIPP::getTravelTime(int start, int end, const ConstraintTable& constraint_table,
                        int upper_bound) {
  reset();
  min_f_val = -1;  // this disables focal list
  int length = MAX_TIMESTEP;
  auto root =
      new SIPPNode(State(start, 0, 0), 0, compute_heuristic(start, end), nullptr, 1, 1, 0, 0);
  pushNodeToOpenAndFocal(root);
  auto static_timestep =
      constraint_table.getMaxTimestep();  // everything is static after this timestep
  while (!open_list.empty()) {
    auto curr = open_list.top();
    open_list.pop();
    if (curr->state.location == end) {
      length = curr->g_val;
      break;
    }
    list<State> next_states = instance.getNeighbors(curr->state);
    next_states.emplace_back(curr->state);
    for (auto& next_state : next_states) {
      int next_timestep = curr->state.timestep + 1;
      int next_g_val = curr->g_val + 1;
      if (static_timestep <= curr->state.timestep) {
        if (curr->state == next_state) {
          continue;
        }
        next_timestep--;
      }
      if (!constraint_table.constrained(next_state.location, next_timestep) &&
          !constraint_table.constrained(curr->state.location, next_state.location, next_timestep)) {
        // if that grid is not blocked
        int next_h_val = compute_heuristic(next_state.location, end);
        if (next_g_val + next_h_val >=
            upper_bound)  // the cost of the path is larger than the upper bound
          continue;
        auto next =
            new SIPPNode(State(next_state.location, next_g_val, next_state.orientation), next_g_val,
                         next_h_val, nullptr, next_timestep + 1, next_timestep + 1, 0, 0);
        if (dominanceCheck(next))
          pushNodeToOpenAndFocal(next);
        else
          delete next;
      }
    }
  }
  releaseNodes();
  num_expanded = 0;
  num_generated = 0;
  num_reopened = 0;
  return length;
}

void SIPP::updateFocalList() {
  auto open_head = open_list.top();
  if (open_head->getFVal() > min_f_val) {
    int new_min_f_val = (int)open_head->getFVal();
    for (auto n : open_list) {
      if (n->getFVal() > w * min_f_val && n->getFVal() <= w * new_min_f_val)
        n->focal_handle = focal_list.push(n);
    }
    min_f_val = new_min_f_val;
  }
}

inline void SIPP::pushNodeToOpenAndFocal(SIPPNode* node) {
  num_generated++;
  node->open_handle = open_list.push(node);
  node->in_openlist = true;
  if (node->getFVal() <= w * min_f_val) node->focal_handle = focal_list.push(node);
  allNodes_table[node].push_back(node);
}

inline void SIPP::pushNodeToFocal(SIPPNode* node) {
  num_generated++;
  allNodes_table[node].push_back(node);
  node->in_openlist = true;
  node->focal_handle = focal_list.push(node);  // we only use focal list; no open list is used
}

inline void SIPP::eraseNodeFromLists(SIPPNode* node) {
  if (open_list.empty()) {
    // we only have focal list
    focal_list.erase(node->focal_handle);
  } else if (focal_list.empty()) {
    // we only have open list
    open_list.erase(node->open_handle);
  } else {
    // we have both open and focal
    open_list.erase(node->open_handle);
    if (node->getFVal() <= w * min_f_val) focal_list.erase(node->focal_handle);
  }
}

void SIPP::releaseNodes() {
  open_list.clear();
  focal_list.clear();
  for (auto& node_list : allNodes_table)
    for (auto n : node_list.second) delete n;
  allNodes_table.clear();
  for (auto n : useless_nodes) delete n;
  useless_nodes.clear();
}

void SIPP::printSearchTree() const {
  vector<list<SIPPNode*>> nodes;
  for (const auto& node_list : allNodes_table) {
    for (const auto& n : node_list.second) {
      if (nodes.size() <= n->state.timestep) nodes.resize(n->state.timestep + 1);
      nodes[n->state.timestep].emplace_back(n);
    }
  }
  cout << "Search Tree" << endl;
  for (int t = 0; t < nodes.size(); t++) {
    cout << "t=" << t << ":\t";
    for (const auto& n : nodes[t])
      cout << *n << "[" << n->state.timestep << "," << n->high_expansion
           << "),c=" << n->num_of_conflicts << "\t";
    cout << endl;
  }
}

// return true iff we the new node is not dominated by any old node
bool SIPP::dominanceCheck(SIPPNode* new_node) {
  auto ptr = allNodes_table.find(new_node);
  if (ptr == allNodes_table.end()) return true;
  for (auto& old_node : ptr->second) {
    if (old_node->state.timestep <= new_node->state.timestep and
        old_node->num_of_conflicts <= new_node->num_of_conflicts) {
      // the new node is dominated by the old node
      return false;
    } else if (old_node->state.timestep >= new_node->state.timestep and
               old_node->num_of_conflicts >=
                   new_node->num_of_conflicts)  // the old node is dominated by the new node
    {
      // delete the old node
      if (old_node->in_openlist)       // the old node has not been expanded yet
        eraseNodeFromLists(old_node);  // delete it from open and/or focal lists
      else                             // the old node has been expanded already
        num_reopened++;                // re-expand it
      useless_nodes.push_back(old_node);
      ptr->second.remove(old_node);
      num_generated--;
      // this is because we later will increase num_generated when we insert the new node into
      // lists.
      return true;
    } else if (old_node->state.timestep < new_node->high_expansion and
               new_node->state.timestep < old_node->high_expansion) {
      // intervals overlap --> we need to split the node to make them disjoint
      if (old_node->state.timestep <= new_node->state.timestep) {
        assert(old_node->num_of_conflicts > new_node->num_of_conflicts);
        old_node->high_expansion = new_node->state.timestep;
      } else  // i.e., old_node->timestep > new_node->timestep
      {
        assert(old_node->num_of_conflicts <= new_node->num_of_conflicts);
        new_node->high_expansion = old_node->state.timestep;
      }
    }
  }
  return true;
}
