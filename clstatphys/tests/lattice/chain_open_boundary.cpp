#include <gtest/gtest.h>
#include <lattice/chain_open_boundary.hpp>

TEST(ChainOpenBoundaryTest, FundamentalTest){
  int Ns = 10;
  lattice::ChainOpenBoundary lattice(Ns);
  int num_particles = lattice.set_num_particles(Ns); 
  ASSERT_EQ(num_particles, Ns);
  int N_adj = lattice.number_adjacent();
  ASSERT_EQ(N_adj,2);
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));;
  lattice.create_table(pair_table);
  int target_site_1 = lattice.numberize(0);
  int target_site_2 = lattice.numberize(num_particles-1);
  EXPECT_EQ(pair_table[target_site_1][0], lattice.numberize(0));
  EXPECT_EQ(pair_table[target_site_1][1], lattice.numberize(1));
  EXPECT_EQ(pair_table[target_site_2][0], lattice.numberize(num_particles-2));
  EXPECT_EQ(pair_table[target_site_2][1], lattice.numberize(num_particles-1));
}


