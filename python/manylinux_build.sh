declare -a arr=("/opt/python/cp36-cp36m/bin" "/opt/python/cp37-cp37m/bin" "/opt/python/cp38-cp38m/bin")

# Compile wheels
for PYBIN in "${arr[@]}"; do
    "${PYBIN}/pip" install -r /io/python/requirements.txt
    "${PYBIN}/pip" wheel /io/python/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in /io/python/wheelhouse/*.whl; do
    auditwheel repair "$whl" -w /io/python/wheelhouse/repaired
done
