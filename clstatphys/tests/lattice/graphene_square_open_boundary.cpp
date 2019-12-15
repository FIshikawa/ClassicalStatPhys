#include <gtest/gtest.h>
#include <lattice/graphene_square_open_boundary.hpp>

TEST(GrapheneSquareOpenBoundaryTest, FundamentalTest){
  int Ns = 3;
  lattice::GrapheneSquareOpenBoundary lattice(Ns);
  int num_particles = lattice.set_num_particles(Ns); 
  ASSERT_EQ(num_particles, 2*Ns*(2*Ns+1));
  int N_adj = lattice.number_adjacent();
  ASSERT_EQ(N_adj,3);
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));;
  lattice.create_table(pair_table);
  ASSERT_EQ(lattice.numberize(0,0),0);
  ASSERT_EQ(lattice.numberize(2*Ns,2*Ns-1),num_particles-1);
  ASSERT_EQ(lattice.latticize(0,num_particles-1),2*Ns);
  ASSERT_EQ(lattice.latticize(1,num_particles-1),2*Ns-1);
  ASSERT_EQ(lattice.target_number(1,1),2);
  int target_site_1 = lattice.numberize(0,0);
  int target_site_2 = lattice.numberize(0,2*Ns-1);
  int target_site_3 = lattice.numberize(2*Ns,2*Ns-1);
  int target_site_4 = lattice.numberize(2*Ns,0);
  int target_site_5 = lattice.numberize(2,2);
  int target_site_6 = lattice.numberize(2,3);
  EXPECT_EQ(pair_table[target_site_1][0], target_site_1);
  EXPECT_EQ(pair_table[target_site_1][1], lattice.numberize(1,0));
  EXPECT_EQ(pair_table[target_site_1][2], lattice.numberize(0,1));
  EXPECT_EQ(pair_table[target_site_2][0], target_site_2);
  EXPECT_EQ(pair_table[target_site_2][1], lattice.numberize(1,2*Ns-1));
  EXPECT_EQ(pair_table[target_site_2][2], lattice.numberize(0,2*Ns-2));
  EXPECT_EQ(pair_table[target_site_3][0], lattice.numberize(2*Ns-1,2*Ns-1));
  EXPECT_EQ(pair_table[target_site_3][1], target_site_3);
  EXPECT_EQ(pair_table[target_site_3][2], lattice.numberize(2*Ns,2*Ns-2));
  EXPECT_EQ(pair_table[target_site_4][0], lattice.numberize(2*Ns-1,0));
  EXPECT_EQ(pair_table[target_site_4][1], target_site_4);
  EXPECT_EQ(pair_table[target_site_4][2], lattice.numberize(2*Ns,1));
  EXPECT_EQ(pair_table[target_site_5][0], lattice.numberize(1,2));
  EXPECT_EQ(pair_table[target_site_5][1], lattice.numberize(3,2));
  EXPECT_EQ(pair_table[target_site_5][2], lattice.numberize(2,3));
  EXPECT_EQ(pair_table[target_site_6][0], lattice.numberize(1,3));
  EXPECT_EQ(pair_table[target_site_6][1], lattice.numberize(3,3));
  EXPECT_EQ(pair_table[target_site_6][2], lattice.numberize(2,2));
}


