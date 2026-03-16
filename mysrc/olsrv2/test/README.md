# RFC5444 Tests

This folder contains a minimal unit test runner for RFC5444 helpers.

## Build and Run

Tests now use mocks for INET dependencies to allow standalone compilation.

### Run NHDP Test
```bash
g++ -std=c++17 -I../include -I. -Imock -DUNIT_TEST \
  test_nhdp.cc \
  ../src/nhdp/nhdp.cc \
  ../src/nhdp/nhdp_db.cc \
  ../src/olsrv2_state.cc \
  -o test_nhdp && ./test_nhdp
```

### Run OLSRv2 Core Test
```bash
g++ -std=c++17 -I../include -I. -Imock -DUNIT_TEST \
  test_olsrv2.cc \
  ../src/olsrv2.cc \
  ../src/olsrv2_routing.cc \
  ../src/olsrv2_tc.cc \
  ../src/olsrv2_state.cc \
  -o test_olsrv2 && ./test_olsrv2
```

### Run Routing Test
```bash
g++ -std=c++17 -I../include -I. -Imock -DUNIT_TEST \
  test_routing.cc \
  ../src/olsrv2_routing.cc \
  ../src/olsrv2_state.cc \
  -o test_routing && ./test_routing
```