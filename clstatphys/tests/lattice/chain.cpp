#include <gtest/gtest.h>
#include <clstatphys/lattice/chain.hpp>

TEST(ChainTest, FundamentalTest){
  int Ns = 10;
  lattice::Chain lattice(Ns);
  int num_particles = lattice.set_num_particles(Ns); 
  ASSERT_EQ(num_particles, Ns);
  int N_adj = lattice.number_adjacent();
  ASSERT_EQ(N_adj,2);
  std::vector<std::vector<int> > pair_table(num_particles,std::vector<int>(N_adj));;
  lattice.create_table(pair_table);
  int target_site = lattice.numberize(0);
  int pair_site_1 = lattice.numberize(num_particles-1);
  int pair_site_2 = lattice.numberize(1);
  EXPECT_EQ(pair_table[target_site][0], pair_site_1);
  EXPECT_EQ(pair_table[target_site][1], pair_site_2);
}


