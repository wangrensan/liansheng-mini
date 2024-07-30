var __getOwnPropNames = Object.getOwnPropertyNames;
var __commonJS = (cb, mod) => function __require() {
  return mod || (0, cb[__getOwnPropNames(cb)[0]])((mod = { exports: {} }).exports, mod), mod.exports;
};

// ../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/detectors/yin.js
var require_yin = __commonJS({
  "../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/detectors/yin.js"(exports) {
    "use strict";
    var __assign = exports && exports.__assign || function() {
      __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
          s = arguments[i];
          for (var p in s)
            if (Object.prototype.hasOwnProperty.call(s, p))
              t[p] = s[p];
        }
        return t;
      };
      return __assign.apply(this, arguments);
    };
    Object.defineProperty(exports, "__esModule", { value: true });
    var DEFAULT_YIN_PARAMS = {
      threshold: 0.1,
      sampleRate: 44100,
      probabilityThreshold: 0.1
    };
    function YIN(params) {
      if (params === void 0) {
        params = {};
      }
      var config = __assign(__assign({}, DEFAULT_YIN_PARAMS), params);
      var threshold = config.threshold, sampleRate = config.sampleRate, probabilityThreshold = config.probabilityThreshold;
      return function YINDetector(float32AudioBuffer) {
        var bufferSize;
        for (bufferSize = 1; bufferSize < float32AudioBuffer.length; bufferSize *= 2)
          ;
        bufferSize /= 2;
        var yinBufferLength = bufferSize / 2;
        var yinBuffer = new Float32Array(yinBufferLength);
        var probability = 0, tau;
        for (var t = 0; t < yinBufferLength; t++) {
          yinBuffer[t] = 0;
        }
        for (var t = 1; t < yinBufferLength; t++) {
          for (var i = 0; i < yinBufferLength; i++) {
            var delta = float32AudioBuffer[i] - float32AudioBuffer[i + t];
            yinBuffer[t] += delta * delta;
          }
        }
        yinBuffer[0] = 1;
        yinBuffer[1] = 1;
        var runningSum = 0;
        for (var t = 1; t < yinBufferLength; t++) {
          runningSum += yinBuffer[t];
          yinBuffer[t] *= t / runningSum;
        }
        for (tau = 2; tau < yinBufferLength; tau++) {
          if (yinBuffer[tau] < threshold) {
            while (tau + 1 < yinBufferLength && yinBuffer[tau + 1] < yinBuffer[tau]) {
              tau++;
            }
            probability = 1 - yinBuffer[tau];
            break;
          }
        }
        if (tau === yinBufferLength || yinBuffer[tau] >= threshold) {
          return null;
        }
        if (probability < probabilityThreshold) {
          return null;
        }
        var betterTau, x0, x2;
        if (tau < 1) {
          x0 = tau;
        } else {
          x0 = tau - 1;
        }
        if (tau + 1 < yinBufferLength) {
          x2 = tau + 1;
        } else {
          x2 = tau;
        }
        if (x0 === tau) {
          if (yinBuffer[tau] <= yinBuffer[x2]) {
            betterTau = tau;
          } else {
            betterTau = x2;
          }
        } else if (x2 === tau) {
          if (yinBuffer[tau] <= yinBuffer[x0]) {
            betterTau = tau;
          } else {
            betterTau = x0;
          }
        } else {
          var s0 = yinBuffer[x0];
          var s1 = yinBuffer[tau];
          var s2 = yinBuffer[x2];
          betterTau = tau + (s2 - s0) / (2 * (2 * s1 - s2 - s0));
        }
        return sampleRate / betterTau;
      };
    }
    exports.YIN = YIN;
  }
});

// ../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/detectors/amdf.js
var require_amdf = __commonJS({
  "../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/detectors/amdf.js"(exports) {
    "use strict";
    var __assign = exports && exports.__assign || function() {
      __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
          s = arguments[i];
          for (var p in s)
            if (Object.prototype.hasOwnProperty.call(s, p))
              t[p] = s[p];
        }
        return t;
      };
      return __assign.apply(this, arguments);
    };
    Object.defineProperty(exports, "__esModule", { value: true });
    var DEFAULT_AMDF_PARAMS = {
      sampleRate: 44100,
      minFrequency: 82,
      maxFrequency: 1e3,
      ratio: 5,
      sensitivity: 0.1
    };
    function AMDF(params) {
      if (params === void 0) {
        params = {};
      }
      var config = __assign(__assign({}, DEFAULT_AMDF_PARAMS), params);
      var sampleRate = config.sampleRate;
      var minFrequency = config.minFrequency;
      var maxFrequency = config.maxFrequency;
      var sensitivity = config.sensitivity;
      var ratio = config.ratio;
      var amd = [];
      var maxPeriod = Math.ceil(sampleRate / minFrequency);
      var minPeriod = Math.floor(sampleRate / maxFrequency);
      return function AMDFDetector(float32AudioBuffer) {
        var maxShift = float32AudioBuffer.length;
        var t = 0;
        var minval = Infinity;
        var maxval = -Infinity;
        var frames1, frames2, calcSub, i, j, u, aux1, aux2;
        for (i = 0; i < maxShift; i++) {
          if (minPeriod <= i && i <= maxPeriod) {
            for (aux1 = 0, aux2 = i, t = 0, frames1 = [], frames2 = []; aux1 < maxShift - i; t++, aux2++, aux1++) {
              frames1[t] = float32AudioBuffer[aux1];
              frames2[t] = float32AudioBuffer[aux2];
            }
            var frameLength = frames1.length;
            calcSub = [];
            for (u = 0; u < frameLength; u++) {
              calcSub[u] = frames1[u] - frames2[u];
            }
            var summation = 0;
            for (u = 0; u < frameLength; u++) {
              summation += Math.abs(calcSub[u]);
            }
            amd[i] = summation;
          }
        }
        for (j = minPeriod; j < maxPeriod; j++) {
          if (amd[j] < minval)
            minval = amd[j];
          if (amd[j] > maxval)
            maxval = amd[j];
        }
        var cutoff = Math.round(sensitivity * (maxval - minval) + minval);
        for (j = minPeriod; j <= maxPeriod && amd[j] > cutoff; j++)
          ;
        var searchLength = minPeriod / 2;
        minval = amd[j];
        var minpos = j;
        for (i = j - 1; i < j + searchLength && i <= maxPeriod; i++) {
          if (amd[i] < minval) {
            minval = amd[i];
            minpos = i;
          }
        }
        if (Math.round(amd[minpos] * ratio) < maxval) {
          return sampleRate / minpos;
        } else {
          return null;
        }
      };
    }
    exports.AMDF = AMDF;
  }
});

// ../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/detectors/acf2plus.js
var require_acf2plus = __commonJS({
  "../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/detectors/acf2plus.js"(exports) {
    "use strict";
    var __assign = exports && exports.__assign || function() {
      __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
          s = arguments[i];
          for (var p in s)
            if (Object.prototype.hasOwnProperty.call(s, p))
              t[p] = s[p];
        }
        return t;
      };
      return __assign.apply(this, arguments);
    };
    Object.defineProperty(exports, "__esModule", { value: true });
    var DEFAULT_PARAMS = {
      sampleRate: 44100
    };
    function ACF2PLUS(params) {
      if (params === void 0) {
        params = DEFAULT_PARAMS;
      }
      var config = __assign(__assign({}, DEFAULT_PARAMS), params);
      var sampleRate = config.sampleRate;
      return function ACF2PLUSDetector(float32AudioBuffer) {
        var maxShift = float32AudioBuffer.length;
        var rms = 0;
        var i, j, u, tmp;
        for (i = 0; i < maxShift; i++) {
          tmp = float32AudioBuffer[i];
          rms += tmp * tmp;
        }
        rms = Math.sqrt(rms / maxShift);
        if (rms < 0.01)
          return -1;
        var aux1 = 0;
        var aux2 = maxShift - 1;
        var thres = 0.2;
        for (i = 0; i < maxShift / 2; i++)
          if (Math.abs(float32AudioBuffer[i]) < thres) {
            aux1 = i;
            break;
          }
        for (i = 1; i < maxShift / 2; i++)
          if (Math.abs(float32AudioBuffer[maxShift - i]) < thres) {
            aux2 = maxShift - i;
            break;
          }
        var frames = float32AudioBuffer.slice(aux1, aux2);
        var framesLength = frames.length;
        var calcSub = new Array(framesLength).fill(0);
        for (i = 0; i < framesLength; i++)
          for (j = 0; j < framesLength - i; j++)
            calcSub[i] = calcSub[i] + frames[j] * frames[j + i];
        u = 0;
        while (calcSub[u] > calcSub[u + 1])
          u++;
        var maxval = -1, maxpos = -1;
        for (i = u; i < framesLength; i++) {
          if (calcSub[i] > maxval) {
            maxval = calcSub[i];
            maxpos = i;
          }
        }
        var T0 = maxpos;
        var y1 = calcSub[T0 - 1], y2 = calcSub[T0], y3 = calcSub[T0 + 1];
        var a = (y1 + y3 - 2 * y2) / 2;
        var b = (y3 - y1) / 2;
        if (a)
          T0 = T0 - b / (2 * a);
        return sampleRate / T0;
      };
    }
    exports.ACF2PLUS = ACF2PLUS;
  }
});

// ../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/detectors/dynamic_wavelet.js
var require_dynamic_wavelet = __commonJS({
  "../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/detectors/dynamic_wavelet.js"(exports) {
    "use strict";
    var __assign = exports && exports.__assign || function() {
      __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
          s = arguments[i];
          for (var p in s)
            if (Object.prototype.hasOwnProperty.call(s, p))
              t[p] = s[p];
        }
        return t;
      };
      return __assign.apply(this, arguments);
    };
    Object.defineProperty(exports, "__esModule", { value: true });
    var MAX_FLWT_LEVELS = 6;
    var MAX_F = 3e3;
    var DIFFERENCE_LEVELS_N = 3;
    var MAXIMA_THRESHOLD_RATIO = 0.75;
    var DEFAULT_DYNAMIC_WAVELET_CONFIG = {
      sampleRate: 44100
    };
    function DynamicWavelet(params) {
      if (params === void 0) {
        params = {};
      }
      var config = __assign(__assign({}, DEFAULT_DYNAMIC_WAVELET_CONFIG), params);
      var sampleRate = config.sampleRate;
      return function DynamicWaveletDetector(float32AudioBuffer) {
        var mins = [];
        var maxs = [];
        var bufferLength = float32AudioBuffer.length;
        var freq = null;
        var theDC = 0;
        var minValue = 0;
        var maxValue = 0;
        for (var i = 0; i < bufferLength; i++) {
          var sample = float32AudioBuffer[i];
          theDC = theDC + sample;
          maxValue = Math.max(maxValue, sample);
          minValue = Math.min(minValue, sample);
        }
        theDC /= bufferLength;
        minValue -= theDC;
        maxValue -= theDC;
        var amplitudeMax = maxValue > -1 * minValue ? maxValue : -1 * minValue;
        var amplitudeThreshold = amplitudeMax * MAXIMA_THRESHOLD_RATIO;
        var curLevel = 0;
        var curModeDistance = -1;
        var curSamNb = float32AudioBuffer.length;
        var delta, nbMaxs, nbMins;
        while (true) {
          delta = ~~(sampleRate / (Math.pow(2, curLevel) * MAX_F));
          if (curSamNb < 2)
            break;
          var dv = void 0;
          var previousDV = -1e3;
          var lastMinIndex = -1e6;
          var lastMaxIndex = -1e6;
          var findMax = false;
          var findMin = false;
          nbMins = 0;
          nbMaxs = 0;
          for (var i = 2; i < curSamNb; i++) {
            var si = float32AudioBuffer[i] - theDC;
            var si1 = float32AudioBuffer[i - 1] - theDC;
            if (si1 <= 0 && si > 0)
              findMax = true;
            if (si1 >= 0 && si < 0)
              findMin = true;
            dv = si - si1;
            if (previousDV > -1e3) {
              if (findMin && previousDV < 0 && dv >= 0) {
                if (Math.abs(si) >= amplitudeThreshold) {
                  if (i > lastMinIndex + delta) {
                    mins[nbMins++] = i;
                    lastMinIndex = i;
                    findMin = false;
                  }
                }
              }
              if (findMax && previousDV > 0 && dv <= 0) {
                if (Math.abs(si) >= amplitudeThreshold) {
                  if (i > lastMaxIndex + delta) {
                    maxs[nbMaxs++] = i;
                    lastMaxIndex = i;
                    findMax = false;
                  }
                }
              }
            }
            previousDV = dv;
          }
          if (nbMins === 0 && nbMaxs === 0) {
            break;
          }
          var d = void 0;
          var distances = [];
          for (var i = 0; i < curSamNb; i++) {
            distances[i] = 0;
          }
          for (var i = 0; i < nbMins; i++) {
            for (var j = 1; j < DIFFERENCE_LEVELS_N; j++) {
              if (i + j < nbMins) {
                d = Math.abs(mins[i] - mins[i + j]);
                distances[d] += 1;
              }
            }
          }
          var bestDistance = -1;
          var bestValue = -1;
          for (var i = 0; i < curSamNb; i++) {
            var summed = 0;
            for (var j = -1 * delta; j <= delta; j++) {
              if (i + j >= 0 && i + j < curSamNb) {
                summed += distances[i + j];
              }
            }
            if (summed === bestValue) {
              if (i === 2 * bestDistance) {
                bestDistance = i;
              }
            } else if (summed > bestValue) {
              bestValue = summed;
              bestDistance = i;
            }
          }
          var distAvg = 0;
          var nbDists = 0;
          for (var j = -delta; j <= delta; j++) {
            if (bestDistance + j >= 0 && bestDistance + j < bufferLength) {
              var nbDist = distances[bestDistance + j];
              if (nbDist > 0) {
                nbDists += nbDist;
                distAvg += (bestDistance + j) * nbDist;
              }
            }
          }
          distAvg /= nbDists;
          if (curModeDistance > -1) {
            if (Math.abs(distAvg * 2 - curModeDistance) <= 2 * delta) {
              freq = sampleRate / (Math.pow(2, curLevel - 1) * curModeDistance);
              break;
            }
          }
          curModeDistance = distAvg;
          curLevel++;
          if (curLevel >= MAX_FLWT_LEVELS || curSamNb < 2) {
            break;
          }
          var newFloat32AudioBuffer = float32AudioBuffer.subarray(0);
          if (curSamNb === distances.length) {
            newFloat32AudioBuffer = new Float32Array(curSamNb / 2);
          }
          for (var i = 0; i < curSamNb / 2; i++) {
            newFloat32AudioBuffer[i] = (float32AudioBuffer[2 * i] + float32AudioBuffer[2 * i + 1]) / 2;
          }
          float32AudioBuffer = newFloat32AudioBuffer;
          curSamNb /= 2;
        }
        return freq;
      };
    }
    exports.DynamicWavelet = DynamicWavelet;
  }
});

// ../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/detectors/macleod.js
var require_macleod = __commonJS({
  "../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/detectors/macleod.js"(exports) {
    "use strict";
    var __assign = exports && exports.__assign || function() {
      __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
          s = arguments[i];
          for (var p in s)
            if (Object.prototype.hasOwnProperty.call(s, p))
              t[p] = s[p];
        }
        return t;
      };
      return __assign.apply(this, arguments);
    };
    Object.defineProperty(exports, "__esModule", { value: true });
    var DEFAULT_MACLEOD_PARAMS = {
      bufferSize: 1024,
      cutoff: 0.97,
      sampleRate: 44100
    };
    function Macleod(params) {
      if (params === void 0) {
        params = {};
      }
      var config = __assign(__assign({}, DEFAULT_MACLEOD_PARAMS), params);
      var bufferSize = config.bufferSize, cutoff = config.cutoff, sampleRate = config.sampleRate;
      var SMALL_CUTOFF = 0.5;
      var LOWER_PITCH_CUTOFF = 80;
      var nsdf = new Float32Array(bufferSize);
      var squaredBufferSum = new Float32Array(bufferSize);
      var turningPointX;
      var turningPointY;
      var maxPositions = [];
      var periodEstimates = [];
      var ampEstimates = [];
      function normalizedSquareDifference(float32AudioBuffer) {
        var acf;
        var divisorM;
        squaredBufferSum[0] = float32AudioBuffer[0] * float32AudioBuffer[0];
        for (var i = 1; i < float32AudioBuffer.length; i += 1) {
          squaredBufferSum[i] = float32AudioBuffer[i] * float32AudioBuffer[i] + squaredBufferSum[i - 1];
        }
        for (var tau = 0; tau < float32AudioBuffer.length; tau++) {
          acf = 0;
          divisorM = squaredBufferSum[float32AudioBuffer.length - 1 - tau] + squaredBufferSum[float32AudioBuffer.length - 1] - squaredBufferSum[tau];
          for (var i = 0; i < float32AudioBuffer.length - tau; i++) {
            acf += float32AudioBuffer[i] * float32AudioBuffer[i + tau];
          }
          nsdf[tau] = 2 * acf / divisorM;
        }
      }
      function parabolicInterpolation(tau) {
        var nsdfa = nsdf[tau - 1], nsdfb = nsdf[tau], nsdfc = nsdf[tau + 1], bValue = tau, bottom = nsdfc + nsdfa - 2 * nsdfb;
        if (bottom === 0) {
          turningPointX = bValue;
          turningPointY = nsdfb;
        } else {
          var delta = nsdfa - nsdfc;
          turningPointX = bValue + delta / (2 * bottom);
          turningPointY = nsdfb - delta * delta / (8 * bottom);
        }
      }
      function peakPicking() {
        var pos = 0;
        var curMaxPos = 0;
        while (pos < (nsdf.length - 1) / 3 && nsdf[pos] > 0) {
          pos++;
        }
        while (pos < nsdf.length - 1 && nsdf[pos] <= 0) {
          pos++;
        }
        if (pos == 0) {
          pos = 1;
        }
        while (pos < nsdf.length - 1) {
          if (nsdf[pos] > nsdf[pos - 1] && nsdf[pos] >= nsdf[pos + 1]) {
            if (curMaxPos == 0) {
              curMaxPos = pos;
            } else if (nsdf[pos] > nsdf[curMaxPos]) {
              curMaxPos = pos;
            }
          }
          pos++;
          if (pos < nsdf.length - 1 && nsdf[pos] <= 0) {
            if (curMaxPos > 0) {
              maxPositions.push(curMaxPos);
              curMaxPos = 0;
            }
            while (pos < nsdf.length - 1 && nsdf[pos] <= 0) {
              pos++;
            }
          }
        }
        if (curMaxPos > 0) {
          maxPositions.push(curMaxPos);
        }
      }
      return function Macleod2(float32AudioBuffer) {
        var pitch;
        maxPositions = [];
        periodEstimates = [];
        ampEstimates = [];
        normalizedSquareDifference(float32AudioBuffer);
        peakPicking();
        var highestAmplitude = -Infinity;
        for (var i = 0; i < maxPositions.length; i++) {
          var tau = maxPositions[i];
          highestAmplitude = Math.max(highestAmplitude, nsdf[tau]);
          if (nsdf[tau] > SMALL_CUTOFF) {
            parabolicInterpolation(tau);
            ampEstimates.push(turningPointY);
            periodEstimates.push(turningPointX);
            highestAmplitude = Math.max(highestAmplitude, turningPointY);
          }
        }
        if (periodEstimates.length) {
          var actualCutoff = cutoff * highestAmplitude;
          var periodIndex = 0;
          for (var i = 0; i < ampEstimates.length; i++) {
            if (ampEstimates[i] >= actualCutoff) {
              periodIndex = i;
              break;
            }
          }
          var period = periodEstimates[periodIndex], pitchEstimate = sampleRate / period;
          if (pitchEstimate > LOWER_PITCH_CUTOFF) {
            pitch = pitchEstimate;
          } else {
            pitch = -1;
          }
        } else {
          pitch = -1;
        }
        return {
          probability: highestAmplitude,
          freq: pitch
        };
      };
    }
    exports.Macleod = Macleod;
  }
});

// ../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/tools/frequencies.js
var require_frequencies = __commonJS({
  "../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/tools/frequencies.js"(exports) {
    "use strict";
    var __assign = exports && exports.__assign || function() {
      __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
          s = arguments[i];
          for (var p in s)
            if (Object.prototype.hasOwnProperty.call(s, p))
              t[p] = s[p];
        }
        return t;
      };
      return __assign.apply(this, arguments);
    };
    Object.defineProperty(exports, "__esModule", { value: true });
    exports.DEFAULT_FREQUENCIES_PARAMS = {
      tempo: 120,
      quantization: 4,
      sampleRate: 44100
    };
    function pitchConsensus(detectors, chunk) {
      var pitches = detectors.map(function(fn) {
        return fn(chunk);
      }).filter(function(value) {
        return value !== null;
      }).sort(function(a, b) {
        return a - b;
      });
      if (pitches.length === 1) {
        return pitches[0];
      } else if (pitches.length === 2) {
        var first = pitches[0], second = pitches[1];
        return first * 2 > second ? Math.sqrt(first * second) : first;
      } else {
        var first = pitches[0];
        var second = pitches[1];
        var secondToLast = pitches[pitches.length - 2];
        var last = pitches[pitches.length - 1];
        var filtered1 = first * 2 > second ? pitches : pitches.slice(1);
        var filtered2 = secondToLast * 2 > last ? filtered1 : filtered1.slice(0, -1);
        return Math.pow(filtered2.reduce(function(t, p) {
          return t * p;
        }, 1), 1 / filtered2.length);
      }
    }
    function frequencies(detector, float32AudioBuffer, options) {
      if (options === void 0) {
        options = {};
      }
      var config = __assign(__assign({}, exports.DEFAULT_FREQUENCIES_PARAMS), options);
      var tempo = config.tempo, quantization = config.quantization, sampleRate = config.sampleRate;
      var bufferLength = float32AudioBuffer.length;
      var chunkSize = Math.round(sampleRate * 60 / (quantization * tempo));
      var getPitch;
      if (Array.isArray(detector)) {
        getPitch = pitchConsensus.bind(null, detector);
      } else {
        getPitch = detector;
      }
      var pitches = [];
      for (var i = 0, max = bufferLength - chunkSize; i <= max; i += chunkSize) {
        var chunk = float32AudioBuffer.slice(i, i + chunkSize);
        var pitch = getPitch(chunk);
        pitches.push(pitch);
      }
      return pitches;
    }
    exports.frequencies = frequencies;
  }
});

// ../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/index.js
var require_lib = __commonJS({
  "../../../Users/18818/Documents/HBuilderProjects/liansheng-mini/node_modules/pitchfinder/lib/index.js"(exports) {
    Object.defineProperty(exports, "__esModule", { value: true });
    var yin_1 = require_yin();
    exports.YIN = yin_1.YIN;
    var amdf_1 = require_amdf();
    exports.AMDF = amdf_1.AMDF;
    var acf2plus_1 = require_acf2plus();
    exports.ACF2PLUS = acf2plus_1.ACF2PLUS;
    var dynamic_wavelet_1 = require_dynamic_wavelet();
    exports.DynamicWavelet = dynamic_wavelet_1.DynamicWavelet;
    var macleod_1 = require_macleod();
    exports.Macleod = macleod_1.Macleod;
    var frequencies_1 = require_frequencies();
    exports.frequencies = frequencies_1.frequencies;
    exports.default = {
      YIN: yin_1.YIN,
      AMDF: amdf_1.AMDF,
      ACF2PLUS: acf2plus_1.ACF2PLUS,
      DynamicWavelet: dynamic_wavelet_1.DynamicWavelet,
      Macleod: macleod_1.Macleod,
      frequencies: frequencies_1.frequencies
    };
  }
});
export default require_lib();
//# sourceMappingURL=pitchfinder.js.map
