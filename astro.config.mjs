// @ts-check
import mdx from '@astrojs/mdx';
import sitemap from '@astrojs/sitemap';
import { defineConfig } from 'astro/config';
import expressiveCode from 'astro-expressive-code';


export default defineConfig({
  image: {
    service: {
      entrypoint: 'astro/assets/services/sharp',
      config: {
        limitInputPixels: false,
      }
    },
    // 设置默认质量
    domains: [],
    remotePatterns: [],
  },
    site: 'https://example.com',
    integrations: [
        expressiveCode({
            // theme
            themes: ['dracula'], 

            styleOverrides: {
                borderRadius: '4px',
                codeFontSize: '1.0rem',
            }
        }),
        mdx(),
        sitemap()
    ],
    markdown: {
        remarkPlugins: [],
        rehypePlugins: [],
    },
    vite: {
        ssr: {
            noExternal: ['marked']
        }
    }
});