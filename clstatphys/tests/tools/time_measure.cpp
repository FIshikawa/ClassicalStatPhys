#include <gtest/gtest.h>
#include <iostream>
#include <tools/time_measure.hpp>
#include<limits>

TEST(TimeMeasureTest, LinearTest){
  double t = 10;
  int N_time = 100;
  double dt = t / N_time;
  int N_time_measure = 10;
  double precision = 1e-6;     

  tools::TimeMeasure time_measure(t, dt, N_time_measure);
  std::vector<double> answer{0.0,1,2,3,4,5,6,7,8,9};
  ASSERT_EQ(time_measure.number_of_total_data(),N_time_measure);
  for(int i = 0; i < time_measure.number_of_total_data(); ++i) 
    std::cout << "i : " << i << ", time_measure[i] : " << time_measure() <<", answer[i] : " << answer[i] << std::endl;
  time_measure.reset();
  double now_time = 0;
  int counter = 0;
  for(int step = 0; step < N_time; ++step){
    if(time_measure.check(step)){
      ASSERT_NEAR(now_time,answer[counter],precision);
      ++counter;
    }
    now_time += dt;
  }
  ASSERT_EQ(counter,time_measure.number_of_total_data());
}

TEST(TimeMeasureTest, LogTest1){
  double t = 10;
  int N_time = 1000;
  double dt = t / N_time;
  int N_time_measure = 10;
  double precision = 1e-6;     

  tools::TimeMeasure time_measure(t, dt, N_time_measure,"log",2);
  std::vector<double> answer{0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9};
  EXPECT_EQ(answer.size(),time_measure.number_of_total_data());
  for(int i = 0; i < time_measure.number_of_total_data(); ++i) 
    std::cout << "i : " << i << ", time_measure[i] : " << time_measure() <<", answer[i] : " << answer[i] << std::endl;
  time_measure.reset();
  double now_time = 0;
  int counter = 0;
  for(int step = 0; step < N_time; ++step){
    if(time_measure.check(step)){
      ASSERT_NEAR(now_time,answer[counter],precision);
      ++counter;
    }
    now_time += dt;
  }
  ASSERT_EQ(counter,time_measure.number_of_total_data());
}

TEST(TimeMeasureTest, LogTest2){
  double t = 20;
  int N_time = 1000;
  double dt = t / N_time;
  int N_time_measure = 10;
  double precision = 1e-6;     

  tools::TimeMeasure time_measure(t, dt, N_time_measure, "log",2);
  std::vector<double> answer{0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2,4,6,8,10,12,14,16,18};
  EXPECT_EQ(answer.size(),time_measure.number_of_total_data());
  for(int i = 0; i < time_measure.number_of_total_data(); ++i) 
    std::cout << "i : " << i << ", time_measure[i] : " << time_measure() <<", answer[i] : " << answer[i] << std::endl;
  time_measure.reset();
  double now_time = 0;
  int counter = 0;
  for(int step = 0; step < N_time; ++step){
    if(time_measure.check(step)){
      ASSERT_NEAR(now_time,answer[counter],precision);
      ++counter;
    }
    now_time += dt;
  }
  ASSERT_EQ(counter,time_measure.number_of_total_data());
}

TEST(TimeMeasureTest, LogTest3){
  double t = 1000;
  int N_time = 100000;
  double dt = t / N_time;
  int N_time_measure = 100;
  double precision = 1e-6;     

  tools::TimeMeasure time_measure(t, dt, N_time_measure, "log");
  std::cout << "order : " << time_measure.order() << std::endl;
  std::cout << "N_total_data : " << time_measure.number_of_total_data() << std::endl;
  for(int i = 0; i < time_measure.number_of_total_data(); ++i) 
    std::cout << "i : " << i << ", time_measure[i] : " << time_measure() << std::endl;
  time_measure.reset();
  double now_time = 0;
  int counter = 0;
  for(int step = 0; step < N_time; ++step){
    if(time_measure.check(step)){
      std::cout <<"counter : " << counter <<", now_time : " << now_time << std::endl; 
      ++counter;
    }
    now_time += dt;
  }
  ASSERT_EQ(counter,time_measure.number_of_total_data());
}
