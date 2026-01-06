# Makefile for Cholesky decomposition project

CC = gcc
CFLAGS = -Wall -Wextra -O2 -std=c11
INCLUDES = -Iinclude
LIBS = -lm

# Directories
SRC_DIR = src
TEST_DIR = tests/cholesky
BUILD_DIR = build

# Source files
CHOLESKY_SRC = $(SRC_DIR)/cholesky.c
TEST_HELPERS_SRC = $(TEST_DIR)/test_helpers.c

# Test executable
TEST_CORRECTNESS = $(BUILD_DIR)/test_correctness

# Default target
all: test-correctness

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Test target
test-correctness: $(BUILD_DIR) $(TEST_CORRECTNESS)
	@echo "Running correctness tests..."
	@$(TEST_CORRECTNESS)

# Build test executable
$(TEST_CORRECTNESS): $(TEST_DIR)/test_correctness.c $(TEST_HELPERS_SRC) $(CHOLESKY_SRC) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(INCLUDES) -I$(TEST_DIR) \
		$(TEST_DIR)/test_correctness.c $(TEST_HELPERS_SRC) $(CHOLESKY_SRC) \
		-o $(TEST_CORRECTNESS) $(LIBS)

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all test-correctness clean

