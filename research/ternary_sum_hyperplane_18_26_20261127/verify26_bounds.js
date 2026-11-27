'use strict';
const assert=require('node:assert/strict');
const H=3n*(3n**13n/339n);
for(const [M,k] of [[3n**15n,6n],[9n*3n**12n,12n],[9n*3n**10n,18n],
    [3n*3n**10n,24n],[H,26n],[H,30n],[3n*1449n,32n]]) {
 assert(M*k**13n<2187n*36n**13n);
}
assert(729+1449<2187);
assert(1449+2*243<2187);
assert(1449+486+243<2187);
assert(2*729>1449);
console.log('PASS: exact all-radius 26D comparisons using published A3(13,6) <= 1449');
