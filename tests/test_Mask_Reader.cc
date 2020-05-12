#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>
#include <sstream>

#include "IBDmix/Mask_Reader.h"

TEST(MaskReader, CanTestInMask) {
  std::istringstream mask_input(
      "1 100 120\n"
      "1 130 140\n"
      "1 160 161\n"
      "1 190 200\n"
      "1 260 281\n"
      "2 130 140\n"
      "4 130 140\n"
      "8 130 140\n"
      "chr8 130 140\n");

  Mask_Reader mask(&mask_input);

  ASSERT_FALSE(mask.in_mask("1", 90));

  // check same position
  ASSERT_FALSE(mask.in_mask("1", 90));
  ASSERT_FALSE(mask.in_mask("1", 90));
  ASSERT_FALSE(mask.in_mask("1", 90));

  ASSERT_FALSE(mask.in_mask("1", 100));

  ASSERT_TRUE(mask.in_mask("1", 101));

  ASSERT_TRUE(mask.in_mask("1", 102));

  ASSERT_TRUE(mask.in_mask("1", 120));

  // skip a range
  ASSERT_TRUE(mask.in_mask("1", 161));

  ASSERT_TRUE(mask.in_mask("1", 261));

  // skip chromosomes
  ASSERT_FALSE(mask.in_mask("2", 130));

  ASSERT_FALSE(mask.in_mask("3", 130));

  ASSERT_TRUE(mask.in_mask("4", 131));

  ASSERT_FALSE(mask.in_mask("5", 131));

  ASSERT_FALSE(mask.in_mask("6", 131));

  ASSERT_TRUE(mask.in_mask("8", 131));

  ASSERT_TRUE(mask.in_mask("chr8", 131));

  ASSERT_FALSE(mask.in_mask("chr9", 131));
  ASSERT_FALSE(mask.in_mask("chr9", 131));
  ASSERT_FALSE(mask.in_mask("chr9", 131));
  ASSERT_FALSE(mask.in_mask("chr9", 131));

  // don't retain the last read values on EOF
  ASSERT_FALSE(mask.in_mask("chr8", 131));

  // setting mask to null should short the region check
  Mask_Reader mask2(nullptr);
  ASSERT_FALSE(mask2.in_mask("1", 161));
}
